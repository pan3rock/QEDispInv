/*
 QEDispInv: Surface Wave Dispersion Curve Computation and Inversion Toolkit
 (MPI Parallelized Version with Optimized Communication)

 GNU General Public License, Version 3, 29 June 2007

 Copyright (c) 2025 Lei Pan

 Xiaofei Chen Research Group,
 Department of Earth and Space Sciences,
 Southern University of Science and Technology, China.
 */

#include "model.hpp"
#include "problem.hpp"
#include "utils.hpp"

#include <CLI11.hpp>
#include <Eigen/Dense>
#include <LBFGSB.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <highfive/H5Easy.hpp>
#include <mpi.h>
#include <random>
#include <string>
#include <toml.hpp>
#include <vector>

using namespace Eigen;
using namespace LBFGSpp;

// Helper function: broadcast a std::string from rank 0 to all other processes
void broadcast_string(std::string &str, int rank) {
  int size;
  if (rank == 0) {
    size = str.length() + 1; // +1 for the null terminator
  }
  // Broadcast the size of the string to all processes
  MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // Create a buffer of the correct size
  std::vector<char> buffer(size);
  if (rank == 0) {
    // If we are rank 0, copy the string data into the buffer
    strcpy(buffer.data(), str.c_str());
  }
  // Broadcast the buffer content to all processes
  MPI_Bcast(buffer.data(), size, MPI_CHAR, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    // If we are not rank 0, create the string from the received buffer
    str = std::string(buffer.data());
  }
}

// Optimization Strategy 1: Serialize all results into a single buffer
std::vector<double> serialize_results(double feval, int it,
                                      const ArrayXd &z_model, const VectorXd &x,
                                      const ArrayXXd &disp) {
  size_t z_size = z_model.size();
  size_t x_size = x.size();
  size_t disp_rows = disp.rows();
  size_t disp_cols = disp.cols();
  size_t disp_size = disp_rows * disp_cols;

  // Buffer layout:
  // [feval, it, z_size, ...z_data, x_size, ...x_data, disp_rows, disp_cols,
  // ...disp_data]
  size_t total_size = 2 + 1 + z_size + 1 + x_size + 2 + disp_size;
  std::vector<double> buffer;
  buffer.reserve(total_size);

  // Pack scalar values
  buffer.push_back(feval);
  buffer.push_back(static_cast<double>(it));

  // Pack z_model
  buffer.push_back(static_cast<double>(z_size));
  buffer.insert(buffer.end(), z_model.data(), z_model.data() + z_size);

  // Pack x (vs_inv)
  buffer.push_back(static_cast<double>(x_size));
  buffer.insert(buffer.end(), x.data(), x.data() + x_size);

  // Pack disp
  buffer.push_back(static_cast<double>(disp_rows));
  buffer.push_back(static_cast<double>(disp_cols));
  buffer.insert(buffer.end(), disp.data(), disp.data() + disp_size);

  return buffer;
}

// Optimization Strategy 1: Deserialize all results from a single buffer
void deserialize_results(const std::vector<double> &buffer, double &feval,
                         int &it, ArrayXd &z_model, ArrayXd &x,
                         ArrayXXd &disp) {
  size_t current_pos = 0;

  // Unpack scalar values
  feval = buffer[current_pos++];
  it = static_cast<int>(buffer[current_pos++]);

  // Unpack z_model
  size_t z_size = static_cast<size_t>(buffer[current_pos++]);
  z_model = Map<ArrayXd>(const_cast<double *>(&buffer[current_pos]), z_size);
  current_pos += z_size;

  // Unpack x (vs_inv)
  size_t x_size = static_cast<size_t>(buffer[current_pos++]);
  x = Map<ArrayXd>(const_cast<double *>(&buffer[current_pos]), x_size);
  current_pos += x_size;

  // Unpack disp
  size_t disp_rows = static_cast<size_t>(buffer[current_pos++]);
  size_t disp_cols = static_cast<size_t>(buffer[current_pos++]);
  disp = Map<ArrayXXd>(const_cast<double *>(&buffer[current_pos]), disp_rows,
                       disp_cols);
}

int main(int argc, char *argv[]) {
  // =========================================================================
  // 1. MPI Initialization
  // =========================================================================
  MPI_Init(&argc, &argv);
  int world_size; // Total number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank; // Rank of the current process (0, 1, 2, ...)
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Define tags for communication
  const int TAG_TASK = 1;
  const int TAG_RESULT = 2;
  const int TASK_TERMINATE = -1;

  // =========================================================================
  // 2. Parse arguments and load data (all processes execute this)
  // =========================================================================
  CLI::App app{"Run the inversion."};
  std::string file_config = "config.toml";
  app.add_option("-c,--config", file_config, "toml-type configure file");
  std::string file_data;
  app.add_option("-d,--data", file_data, "filename of dispersion curves");
  std::string file_mref;
  app.add_option("--model_ref", file_mref, "filename of reference model");
  bool sh = false;
  app.add_flag("--sh", sh, "whether are Love waves");
  std::string file_out = "inv.h5";
  app.add_option("-o,--out", file_out, "filename of output");

  if (world_rank == 0) {
    try {
      app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
      int result = app.exit(e);
      MPI_Abort(MPI_COMM_WORLD, result);
    }
  }

  broadcast_string(file_config, world_rank);
  broadcast_string(file_data, world_rank);
  broadcast_string(file_mref, world_rank);
  broadcast_string(file_out, world_rank);
  MPI_Bcast(&sh, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

  const auto config = toml::parse(file_config);
  auto conf_inv = toml::find(config, "inversion");

  if (file_mref == "")
    file_mref = toml::find<std::string>(conf_inv, "model_ref");
  const auto mref_type =
      toml::find_or<std::string>(conf_inv, "mref_type", "linear");
  auto vs2model = toml::find<std::string>(conf_inv, "vs2model");
  const auto vs_width = toml::find<double>(conf_inv, "vs_width");
  const auto lamb_vs = toml::find<double>(conf_inv, "lambda");
  const auto rtype = toml::find<int>(conf_inv, "reg_type");
  const auto weight = toml::find<std::vector<double>>(conf_inv, "weight");

  ArrayXXd model_ref = loadtxt(file_mref);
  ArrayXXd data_input = loadtxt(file_data);

  std::transform(vs2model.begin(), vs2model.end(), vs2model.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  std::shared_ptr<Vs2Model> pmodel;
  if (vs2model == "nearsurface") {
    pmodel = std::make_shared<NearSurface>(model_ref, mref_type);
    const auto vp2vs = toml::find<double>(conf_inv, "vp2vs");
    pmodel->set_param({vp2vs});
  } else if (vs2model == "gardner") {
    pmodel = std::make_shared<Gardner>(model_ref, mref_type);
  } else if (vs2model == "fixvprho") {
    pmodel = std::make_shared<FixVpRho>(model_ref, mref_type);
  } else if (vs2model == "brocher05") {
    pmodel = std::make_shared<Brocher05>(model_ref, mref_type);
  } else {
    if (world_rank == 0) {
      std::string msg =
          fmt::format("invalid vs2model: {:s}, not in (NearSurface, Gardner, "
                      "FixVpRho, Brocher05)",
                      vs2model);
      std::cerr << msg << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  LBFGSBParam<double> param;
  param.max_iterations = 100;
  param.epsilon = 1e-4;
  param.past = 5;
  param.delta = 1e-6;
  LBFGSBSolver<double> solver(param);

  const auto num_init = toml::find<int>(conf_inv, "num_init");
  const auto num_noise = toml::find<int>(conf_inv, "num_noise");
  const auto rand_depth = toml::find<bool>(conf_inv, "rand_depth");
  const auto rand_vs = toml::find<bool>(conf_inv, "rand_vs");
  auto zmax = toml::find<double>(conf_inv, "zmax");
  const auto r0 = toml::find<double>(conf_inv, "r0");
  const auto rmin = toml::find<double>(conf_inv, "rmin");
  const auto rmax = toml::find<double>(conf_inv, "rmax");
  const auto fix_wavelen = toml::find_or<bool>(conf_inv, "fix_wavelen", false);

  Data data(data_input);
  double lmin, lmax;
  if (fix_wavelen) {
    lmin = toml::find<double>(conf_inv, "lmin");
    lmax = toml::find<double>(conf_inv, "lmax");
  } else {
    lmin = data.lmin;
    lmax = data.lmax;
  }

  if (num_noise > 1 && data_input.cols() < 4) {
    const auto sigma = toml::find<std::vector<double>>(conf_inv, "sigma");
    if (sigma.size() < weight.size()) {
      if (world_rank == 0) {
        std::string msg = fmt::format(
            "the size of sigma ({:d}) is less than that of weight ({:d})",
            sigma.size(), weight.size());
        std::cerr << msg << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    data.add_sigma(sigma);
  }

  std::vector<Data> data_noise;
  if (num_noise == 1) {
    data_noise.push_back(data);
  } else {
    for (int i_d = 0; i_d < num_noise; ++i_d) {
      Data data_resampled = resample(data);
      data_noise.push_back(data_resampled);
    }
  }

  if (world_rank == 0) {
    fmt::print("lmin={:12.3f}, lmax={:12.3f}\n", lmin, lmax);
  }

  std::vector<ArrayXd> z_init;
  if (num_init == 1) {
    z_init.push_back(model_ref.col(1));
  } else {
    for (int i_m = 0; i_m < num_init; ++i_m) {
      if (rand_depth) {
        z_init.push_back(
            generate_depth_by_layer_ratio(lmin, lmax, r0, rmin, rmax, zmax));
      } else {
        z_init.push_back(model_ref.col(1));
      }
    }
  }

  std::vector<ArrayXd> vs_ref, vs_lb, vs_ub;
  for (int i_m = 0; i_m < num_init; ++i_m) {
    ArrayXd z = z_init[i_m];
    int nl = z.rows();
    ArrayXd vsr(nl), lb(nl), ub(nl);
    pmodel->get_vs_limits(z, vs_width, vsr, lb, ub);
    vs_ref.push_back(vsr);
    vs_lb.push_back(lb);
    vs_ub.push_back(ub);
  }

  bool add_last_layer = true;
  if (add_last_layer) {
    auto c0 = data.c[0];
    double cmax0 = *std::max_element(c0.begin(), c0.end());
    double half_lmax = data.lmax / 2.0;
    if (half_lmax > zmax) {
      for (int i_m = 0; i_m < num_init; ++i_m) {
        int nl = z_init[i_m].rows();
        if (half_lmax < z_init[i_m](nl - 1))
          continue;
        double cmax = std::max(cmax0, vs_ref[i_m][nl - 1]);
        ArrayXd tmp(nl + 1);
        tmp.head(nl) = z_init[i_m];
        tmp(nl) = half_lmax;
        z_init[i_m] = tmp;
        tmp.head(nl) = vs_ref[i_m];
        tmp(nl) = cmax * 1.01;
        vs_ref[i_m] = tmp;
        tmp.head(nl) = vs_lb[i_m];
        tmp(nl) = cmax;
        vs_lb[i_m] = tmp;
        tmp.head(nl) = vs_ub[i_m];
        tmp(nl) = cmax + vs_width;
        vs_ub[i_m] = tmp;
      }
    }
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0.0, 1.0);

  auto gen_rand_vsinit = [&](const ArrayXd &lb, const ArrayXd &ub) -> ArrayXd {
    int nl = lb.rows();
    ArrayXd rand = ArrayXd::NullaryExpr(nl, [&]() { return dist(gen); });
    ArrayXd x = lb.array() + (ub - lb).array() * rand;
    return x;
  };

  std::vector<ArrayXd> vs_init;
  if (num_init == 1) {
    vs_init.push_back(vs_ref[0]);
  } else {
    for (int i_m = 0; i_m < num_init; ++i_m) {
      if (rand_vs) {
        vs_init.push_back(gen_rand_vsinit(vs_lb[i_m], vs_ub[i_m]));
      } else {
        auto vs = pmodel->interp_vs(z_init[i_m]);
        vs_init.push_back(vs);
      }
    }
  }

  int num_total = num_init * num_noise;

  // =========================================================================
  // 3. Main parallel logic: Master/Worker Dynamic Scheduling
  // =========================================================================

  if (world_rank == 0) {
    // --- MASTER PROCESS (rank 0) ---
    int next_task_id = 0;
    int tasks_completed = 0;
    MPI_Status status;

    std::vector<ArrayXd> z_inv, vs_inv;
    std::vector<double> fitness;
    std::vector<int> niter;
    std::vector<ArrayXXd> disp_syn;

    z_inv.resize(num_total);
    vs_inv.resize(num_total);
    fitness.resize(num_total);
    niter.resize(num_total);
    disp_syn.resize(num_total);

    // Send initial tasks to all worker processes
    for (int rank = 1; rank < world_size; ++rank) {
      if (next_task_id < num_total) {
        MPI_Send(&next_task_id, 1, MPI_INT, rank, TAG_TASK, MPI_COMM_WORLD);
        next_task_id++;
      } else {
        MPI_Send(&TASK_TERMINATE, 1, MPI_INT, rank, TAG_TASK, MPI_COMM_WORLD);
      }
    }

    // Receive results and distribute new tasks
    while (tasks_completed < num_total) {

      int received_task_id;
      MPI_Recv(&received_task_id, 1, MPI_INT, MPI_ANY_SOURCE, TAG_RESULT,
               MPI_COMM_WORLD, &status);
      int source_rank = status.MPI_SOURCE;

      // Receive the serialized result buffer
      int buffer_size;
      MPI_Recv(&buffer_size, 1, MPI_INT, source_rank, TAG_RESULT,
               MPI_COMM_WORLD, &status);
      std::vector<double> result_buffer(buffer_size);
      MPI_Recv(result_buffer.data(), buffer_size, MPI_DOUBLE, source_rank,
               TAG_RESULT, MPI_COMM_WORLD, &status);

      // Deserialize the buffer and store the results
      double feval;
      int it;
      ArrayXd z_model, x;
      ArrayXXd disp;
      deserialize_results(result_buffer, feval, it, z_model, x, disp);

      fitness[received_task_id] = feval;
      niter[received_task_id] = it;
      z_inv[received_task_id] = z_model;
      vs_inv[received_task_id] = x;
      disp_syn[received_task_id] = disp;

      tasks_completed++;

      // Send the next available task to the worker that just finished
      if (next_task_id < num_total) {
        MPI_Send(&next_task_id, 1, MPI_INT, source_rank, TAG_TASK,
                 MPI_COMM_WORLD);
        next_task_id++;
      } else {
        MPI_Send(&TASK_TERMINATE, 1, MPI_INT, source_rank, TAG_TASK,
                 MPI_COMM_WORLD);
      }
    }

    // =========================================================================
    // 4. Rank 0 performs post-processing and file writing
    // =========================================================================
    if (num_total > 1) {
      std::vector<size_t> idx_outlier = detect_outliers(fitness);
      remove_by_indices(fitness, idx_outlier);
      remove_by_indices(niter, idx_outlier);
      remove_by_indices(z_inv, idx_outlier);
      remove_by_indices(vs_inv, idx_outlier);
      remove_by_indices(disp_syn, idx_outlier);
    }

    const int num_hist = 100;
    double vsmin = min_varray(vs_inv);
    double vsmax = max_varray(vs_inv);
    vsmin *= 0.95;
    vsmax *= 1.05;
    ArrayXd z_samples(num_hist), vs_samples(num_hist);
    ArrayXXd hist = compute_hist2d(z_inv, vs_inv, fitness, vsmin, vsmax, zmax,
                                   num_hist, z_samples, vs_samples);

    ArrayXd vs_mean(num_hist), vs_mode(num_hist), vs_median(num_hist),
        vs_cred10(num_hist), vs_cred90(num_hist);
    compute_statistics(z_samples, vs_samples, hist, vs_mean, vs_median, vs_mode,
                       vs_cred10, vs_cred90);

    H5Easy::File out_h5(file_out, H5Easy::File::Overwrite);
    H5Easy::dump(out_h5, "fitness", fitness);
    H5Easy::dump(out_h5, "niter", niter);
    H5Easy::dump(out_h5, "z_sample", z_samples);
    H5Easy::dump(out_h5, "vs_sample", vs_samples);
    H5Easy::dump(out_h5, "vs_hist2d", hist);
    H5Easy::dump(out_h5, "data", data_input);
    H5Easy::dump(out_h5, "vs_mean", vs_mean);
    H5Easy::dump(out_h5, "vs_median", vs_median);
    H5Easy::dump(out_h5, "vs_mode", vs_mode);
    H5Easy::dump(out_h5, "vs_cred10", vs_cred10);
    H5Easy::dump(out_h5, "vs_cred90", vs_cred90);

    ArrayXXd model_mean = pmodel->generate(z_samples, vs_mean);
    H5Easy::dump(out_h5, "model_mean", model_mean);
    ArrayXd vs_ref_save = pmodel->interp_vs(z_samples);
    H5Easy::dump(out_h5, "vs_ref", vs_ref_save);

    H5Easy::dump(out_h5, "num_init", num_init);
    for (size_t i = 0; i < vs_init.size(); ++i) {
      auto model = pmodel->generate(z_init[i], vs_init[i]);
      std::string key = fmt::format("model_init/{:d}", i);
      H5Easy::dump(out_h5, key, model);
    }

    std::vector<int> mode_used;
    for (size_t i = 0; i < weight.size(); ++i) {
      if (weight[i] > 0) {
        mode_used.push_back(i);
      }
    }
    H5Easy::dump(out_h5, "mode_used", mode_used);

    int num_valid = fitness.size();
    H5Easy::dump(out_h5, "num_valid", num_valid);
    for (int i = 0; i < num_valid; ++i) {
      std::string key = fmt::format("disp/{:d}", i);
      H5Easy::dump(out_h5, key, disp_syn[i]);
    }

  } else {
    // --- WORKER PROCESSES (rank > 0) ---
    MPI_Status status;
    while (true) {
      int task_id;
      MPI_Recv(&task_id, 1, MPI_INT, 0, TAG_TASK, MPI_COMM_WORLD, &status);

      if (task_id == TASK_TERMINATE) {
        break; // Master said no more tasks
      }

      // --- Perform the calculation for the received task_id ---
      int i_d = task_id / num_init;
      int i_m = task_id % num_init;
      auto data_resampled = data_noise[i_d];
      ArrayXd z_model = z_init[i_m];
      lbfgspp::DispersionCurves prob(z_model, vs_ref[i_m], pmodel, weight,
                                     lamb_vs, sh, rtype);
      prob.load_data(data_resampled);

      VectorXd x = vs_init[i_m];
      double feval = 0.0;
      int it = 0;
      try {
        it = solver.minimize(prob, x, feval, vs_lb[i_m], vs_ub[i_m]);
      } catch (const std::exception &exc) {
        // Errors are handled locally and won't crash the program
      }

      // --- Send results back to the master in a single serialized buffer ---
      auto model = pmodel->generate(z_model, x);
      auto disp = prob.forward(model);
      std::vector<double> result_buffer =
          serialize_results(feval, it, z_model, x, disp);
      int buffer_size = result_buffer.size();

      MPI_Send(&task_id, 1, MPI_INT, 0, TAG_RESULT, MPI_COMM_WORLD);
      MPI_Send(&buffer_size, 1, MPI_INT, 0, TAG_RESULT, MPI_COMM_WORLD);
      MPI_Send(result_buffer.data(), buffer_size, MPI_DOUBLE, 0, TAG_RESULT,
               MPI_COMM_WORLD);
    }
  }

  // =========================================================================
  // 5. Finalize MPI
  // =========================================================================
  MPI_Finalize();

  return 0;
}
