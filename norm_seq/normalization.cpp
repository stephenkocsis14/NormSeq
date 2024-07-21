#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <cmath>
#include <vector>
#include <algorithm>

namespace py = pybind11;

py::array_t<double> cpm(py::array_t<int> raw_counts) {
    auto buf = raw_counts.request();
    int n_rows = buf.shape[0];
    int n_cols = buf.shape[1];
    auto ptr = static_cast<int *>(buf.ptr);

    std::vector<double> normalized_counts(n_rows * n_cols);

    for (int j = 0; j < n_cols; ++j) {
        double total_counts = 0;
        for (int i = 0; i < n_rows; ++i) {
            total_counts += ptr[i * n_cols + j];
        }
        for (int i = 0; i < n_rows; ++i) {
            normalized_counts[i * n_cols + j] = (ptr[i * n_cols + j] / total_counts) * 1e6;
        }
    }

    return py::array_t<double>({n_rows, n_cols}, normalized_counts.data());
}

py::array_t<double> tpm(py::array_t<int> raw_counts, py::array_t<double> lengths) {
    auto counts_buf = raw_counts.request();
    auto lengths_buf = lengths.request();
    int n_rows = counts_buf.shape[0];
    int n_cols = counts_buf.shape[1];
    auto counts_ptr = static_cast<int *>(counts_buf.ptr);
    auto lengths_ptr = static_cast<double *>(lengths_buf.ptr);

    std::vector<double> normalized_counts(n_rows * n_cols);
    std::vector<double> rpk(n_rows * n_cols);

    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            rpk[i * n_cols + j] = counts_ptr[i * n_cols + j] / lengths_ptr[i];
        }
    }

    for (int j = 0; j < n_cols; ++j) {
        double sum_rpk = 0;
        for (int i = 0; i < n_rows; ++i) {
            sum_rpk += rpk[i * n_cols + j];
        }
        for (int i = 0; i < n_rows; ++i) {
            normalized_counts[i * n_cols + j] = (rpk[i * n_cols + j] / sum_rpk) * 1e6;
        }
    }

    return py::array_t<double>({n_rows, n_cols}, normalized_counts.data());
}

py::array_t<double> fpkm(py::array_t<int> raw_counts, py::array_t<double> lengths) {
    auto counts_buf = raw_counts.request();
    auto lengths_buf = lengths.request();
    int n_rows = counts_buf.shape[0];
    int n_cols = counts_buf.shape[1];
    auto counts_ptr = static_cast<int *>(counts_buf.ptr);
    auto lengths_ptr = static_cast<double *>(lengths_buf.ptr);

    std::vector<double> normalized_counts(n_rows * n_cols);

    for (int j = 0; j < n_cols; ++j) {
        double total_counts = 0;
        for (int i = 0; i < n_rows; ++i) {
            total_counts += counts_ptr[i * n_cols + j];
        }
        for (int i = 0; i < n_rows; ++i) {
            double per_kilobase = counts_ptr[i * n_cols + j] / (lengths_ptr[i] / 1000.0);
            normalized_counts[i * n_cols + j] = per_kilobase / (total_counts / 1e6);
        }
    }

    return py::array_t<double>({n_rows, n_cols}, normalized_counts.data());
}

py::array_t<double> rpkm(py::array_t<int> raw_counts, py::array_t<double> lengths) {
    // RPKM is essentially the same as FPKM
    return fpkm(raw_counts, lengths);
}

py::array_t<double> log2cpm(py::array_t<int> raw_counts) {
    auto buf = raw_counts.request();
    int n_rows = buf.shape[0];
    int n_cols = buf.shape[1];
    auto ptr = static_cast<int *>(buf.ptr);

    std::vector<double> normalized_counts(n_rows * n_cols);

    for (int j = 0; j < n_cols; ++j) {
        double total_counts = 0;
        for (int i = 0; i < n_rows; ++i) {
            total_counts += ptr[i * n_cols + j];
        }
        for (int i = 0; i < n_rows; ++i) {
            double cpm_value = (ptr[i * n_cols + j] / total_counts) * 1e6;
            normalized_counts[i * n_cols + j] = std::log2(cpm_value + 1);
        }
    }

    return py::array_t<double>({n_rows, n_cols}, normalized_counts.data());
}

py::array_t<double> zscore(py::array_t<int> raw_counts) {
    auto buf = raw_counts.request();
    int n_rows = buf.shape[0];
    int n_cols = buf.shape[1];
    auto ptr = static_cast<int *>(buf.ptr);

    std::vector<double> normalized_counts(n_rows * n_cols);

    for (int j = 0; j < n_cols; ++j) {
        double mean = 0;
        double sq_sum = 0;
        for (int i = 0; i < n_rows; ++i) {
            mean += ptr[i * n_cols + j];
            sq_sum += ptr[i * n_cols + j] * ptr[i * n_cols + j];
        }
        mean /= n_rows;
        double variance = (sq_sum / n_rows) - (mean * mean);
        double stddev = std::sqrt(variance);

        for (int i = 0; i < n_rows; ++i) {
            normalized_counts[i * n_cols + j] = (ptr[i * n_cols + j] - mean) / stddev;
        }
    }

    return py::array_t<double>({n_rows, n_cols}, normalized_counts.data());
}

PYBIND11_MODULE(normalization, m) {
    m.def("cpm", &cpm, "Normalize counts per million (CPM)");
    m.def("tpm", &tpm, "Normalize transcripts per kilobase million (TPM)");
    m.def("fpkm", &fpkm, "Normalize fragments per kilobase million (FPKM)");
    m.def("rpkm", &rpkm, "Normalize reads per kilobase million (RPKM)");
    m.def("log2cpm", &log2cpm, "Log2 normalize counts per million (log2CPM)");
    m.def("zscore", &zscore, "Z-score normalization");
}
