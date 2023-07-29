#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

py::array_t<double> cpp_scaling_2d(py::array_t<double> sigma) {
  py::buffer_info info = sigma.request();
  double* data = static_cast<double*>(info.ptr);
  
  int k = info.shape[1];
  int pz_dim = 1;
  for (int i = 0; i < k; ++i) {
    pz_dim *= 2;
  }
  
  py::array_t<double> pz(pz_dim);
  py::buffer_info pz_info = pz.request();
  double* pz_data = static_cast<double*>(pz_info.ptr);
  
  for (int i = 0; i < pz_dim; ++i) {
    for (int ki = 0; ki < k; ++ki) {
      for (int kj = 0; kj < k; ++kj) {
        pz_data[i] += data[ki * k + kj] * ((i & (1 << ki)) > 0) * ((i & (1 << kj)) > 0);
      }
    }
    pz_data[i] = exp(pz_data[i]);
  }
  
  return pz;
}

PYBIND11_MODULE(cpp_scaling_2d, m) {
  m.def("cpp_scaling_2d", &cpp_scaling_2d, "Compute pz using C++ implementation");
}