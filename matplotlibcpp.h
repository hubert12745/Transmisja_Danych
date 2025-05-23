#pragma once
//
// Wishlist:
// * (WIP) Make work for Eigen Vectors and Matrices
// * export functions in different files for better structure
// * make plot(y) work with x of unsigned type, or get to the bottom of that
//   problem at least
// * errorbar for xerr and yerr
// * errorbar for yerr of shape (2, N)
//
// Changed:
// * Implement a better way for named_plot, maybe just as additional
//   method with extra keyword
// * add location keyword for legend
//

#pragma once

#include <algorithm>
#include <array>
#include <cstdint> // <cstdint> requires c++11 support
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <sstream>

#include <Python.h>

#ifndef WITHOUT_NUMPY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#endif // WITHOUT_NUMPY

#if PY_MAJOR_VERSION >= 3
#define PyString_FromString PyUnicode_FromString
#define PyInt_FromLong PyLong_FromLong
#define PyString_FromString PyUnicode_FromString
#endif

namespace matplotlibcpp {
    namespace detail {

        static std::string s_backend;

        struct _interpreter {
            PyObject* s_python_function_show;
            PyObject* s_python_function_close;
            PyObject* s_python_function_draw;
            PyObject* s_python_function_pause;
            PyObject* s_python_function_save;
            PyObject* s_python_function_figure;
            PyObject* s_python_function_fignum_exists;
            PyObject* s_python_function_plot;
            PyObject* s_python_function_quiver;
            PyObject* s_python_function_contour;
            PyObject* s_python_function_colormap;
            PyObject* s_python_function_axhline;
            PyObject* s_python_function_axvline;
            PyObject* s_python_function_semilogx;
            PyObject* s_python_function_semilogy;
            PyObject* s_python_function_loglog;
            PyObject* s_python_function_fill;
            PyObject* s_python_function_fill_between;
            PyObject* s_python_function_hist;
            PyObject* s_python_function_scatter;
            PyObject* s_python_function_spy;
            PyObject* s_python_function_subplot;
            PyObject* s_python_function_legend;
            PyObject* s_python_function_xlim;
            PyObject* s_python_function_ion;
            PyObject* s_python_function_ginput;
            PyObject* s_python_function_ylim;
            PyObject* s_python_function_title;
            PyObject* s_python_function_axis;
            PyObject* s_python_function_xlabel;
            PyObject* s_python_function_ylabel;
            PyObject* s_python_function_xticks;
            PyObject* s_python_function_yticks;
            PyObject* s_python_function_xscale;
            PyObject* s_python_function_yscale;
            PyObject* s_python_function_grid;
            PyObject* s_python_function_clf;
            PyObject* s_python_function_errorbar;
            PyObject* s_python_function_annotate;
            PyObject* s_python_function_tight_layout;
            PyObject* s_python_colormap;
            PyObject* s_python_empty_tuple;
            PyObject* s_python_function_stem;
            PyObject* s_python_function_xkcd;
            PyObject* s_python_function_text;
            PyObject* s_python_function_suptitle;
            PyObject* s_python_function_bar;
            PyObject* s_python_function_subplots_adjust;
            PyObject* s_python_function_imshow;
            PyObject* s_python_function_colorbar;

            /* For now, _interpreter is implemented as a singleton since its currently not
               possible to have multiple independent embedded python interpreters without
               patching the python source code or starting a separate process for each.
                http://bytes.com/topic/python/answers/793370-multiple-independent-python-interpreters-c-c-program
               */

            static _interpreter& get() {
                static _interpreter ctx;
                return ctx;
            }

        private:
#ifndef WITHOUT_NUMPY
#if PY_MAJOR_VERSION >= 3

            void* import_numpy() {
                import_array(); // initialize C-API
                return NULL;
            }

#else

            void import_numpy() {
                import_array(); // initialize C-API
            }

#endif
#endif

            _interpreter() {

                // optional but recommended
#if PY_MAJOR_VERSION >= 3
                wchar_t name[] = L"plotting";
#else
                char name[] = "plotting";
#endif
                Py_SetProgramName(name);
                Py_Initialize();

#ifndef WITHOUT_NUMPY
                import_numpy(); // initialize numpy C-API
#endif

                PyObject* matplotlibname = PyString_FromString("matplotlib");
                PyObject* pyplotname = PyString_FromString("matplotlib.pyplot");
                PyObject* cmname = PyString_FromString("matplotlib.cm");
                PyObject* pylabname = PyString_FromString("pylab");
                if (!pyplotname || !pylabname || !matplotlibname || !cmname) {
                    throw std::runtime_error("couldnt create string");
                }

                PyObject* matplotlib = PyImport_Import(matplotlibname);
                Py_DECREF(matplotlibname);
                if (!matplotlib) {
                    PyErr_Print();
                    throw std::runtime_error("Error loading module matplotlib!");
                }

                // matplotlib.use() must be called *before* pylab, matplotlib.pyplot,
                // or matplotlib.backends is imported for the first time
                if (!s_backend.empty()) {
                    PyObject_CallMethod(matplotlib, const_cast<char*>("use"),
                        const_cast<char*>("s"), s_backend.c_str());
                }

                PyObject* pymod = PyImport_Import(pyplotname);
                Py_DECREF(pyplotname);
                if (!pymod) {
                    throw std::runtime_error("Error loading module matplotlib.pyplot!");
                }

                s_python_colormap = PyImport_Import(cmname);
                Py_DECREF(cmname);
                if (!s_python_colormap) {
                    throw std::runtime_error("Error loading module matplotlib.cm!");
                }

                PyObject* pylabmod = PyImport_Import(pylabname);
                Py_DECREF(pylabname);
                if (!pylabmod) {
                    throw std::runtime_error("Error loading module pylab!");
                }

                s_python_function_show = PyObject_GetAttrString(pymod, "show");
                s_python_function_close = PyObject_GetAttrString(pymod, "close");
                s_python_function_draw = PyObject_GetAttrString(pymod, "draw");
                s_python_function_pause = PyObject_GetAttrString(pymod, "pause");
                s_python_function_figure = PyObject_GetAttrString(pymod, "figure");
                s_python_function_fignum_exists =
                    PyObject_GetAttrString(pymod, "fignum_exists");
                s_python_function_plot = PyObject_GetAttrString(pymod, "plot");
                s_python_function_quiver = PyObject_GetAttrString(pymod, "quiver");
                s_python_function_contour = PyObject_GetAttrString(pymod, "contour");
                s_python_function_axhline = PyObject_GetAttrString(pymod, "axhline");
                s_python_function_axvline = PyObject_GetAttrString(pymod, "axvline");
                s_python_function_semilogx = PyObject_GetAttrString(pymod, "semilogx");
                s_python_function_semilogy = PyObject_GetAttrString(pymod, "semilogy");
                s_python_function_loglog = PyObject_GetAttrString(pymod, "loglog");
                s_python_function_fill = PyObject_GetAttrString(pymod, "fill");
                s_python_function_fill_between =
                    PyObject_GetAttrString(pymod, "fill_between");
                s_python_function_hist = PyObject_GetAttrString(pymod, "hist");
                s_python_function_scatter = PyObject_GetAttrString(pymod, "scatter");
                s_python_function_spy = PyObject_GetAttrString(pymod, "spy");
                s_python_function_subplot = PyObject_GetAttrString(pymod, "subplot");
                s_python_function_legend = PyObject_GetAttrString(pymod, "legend");
                s_python_function_ylim = PyObject_GetAttrString(pymod, "ylim");
                s_python_function_title = PyObject_GetAttrString(pymod, "title");
                s_python_function_axis = PyObject_GetAttrString(pymod, "axis");
                s_python_function_xlabel = PyObject_GetAttrString(pymod, "xlabel");
                s_python_function_ylabel = PyObject_GetAttrString(pymod, "ylabel");
                s_python_function_xticks = PyObject_GetAttrString(pymod, "xticks");
                s_python_function_yticks = PyObject_GetAttrString(pymod, "yticks");
                s_python_function_xscale = PyObject_GetAttrString(pymod, "xscale");
                s_python_function_yscale = PyObject_GetAttrString(pymod, "yscale");
                s_python_function_grid = PyObject_GetAttrString(pymod, "grid");
                s_python_function_xlim = PyObject_GetAttrString(pymod, "xlim");
                s_python_function_ion = PyObject_GetAttrString(pymod, "ion");
                s_python_function_ginput = PyObject_GetAttrString(pymod, "ginput");
                s_python_function_save = PyObject_GetAttrString(pylabmod, "savefig");
                s_python_function_annotate = PyObject_GetAttrString(pymod, "annotate");
                s_python_function_clf = PyObject_GetAttrString(pymod, "clf");
                s_python_function_errorbar = PyObject_GetAttrString(pymod, "errorbar");
                s_python_function_tight_layout =
                    PyObject_GetAttrString(pymod, "tight_layout");
                s_python_function_stem = PyObject_GetAttrString(pymod, "stem");
                s_python_function_xkcd = PyObject_GetAttrString(pymod, "xkcd");
                s_python_function_text = PyObject_GetAttrString(pymod, "text");
                s_python_function_suptitle = PyObject_GetAttrString(pymod, "suptitle");
                s_python_function_bar = PyObject_GetAttrString(pymod, "bar");
                s_python_function_subplots_adjust =
                    PyObject_GetAttrString(pymod, "subplots_adjust");
                s_python_function_imshow = PyObject_GetAttrString(pymod, "imshow");
                s_python_function_colorbar = PyObject_GetAttrString(pymod, "colorbar");

                if (!s_python_function_show || !s_python_function_close ||
                    !s_python_function_draw || !s_python_function_pause ||
                    !s_python_function_figure || !s_python_function_fignum_exists ||
                    !s_python_function_plot || !s_python_function_quiver ||
                    !s_python_function_contour || !s_python_function_colorbar ||
                    !s_python_function_semilogx || !s_python_function_semilogy ||
                    !s_python_function_loglog || !s_python_function_fill ||
                    !s_python_function_fill_between || !s_python_function_subplot ||
                    !s_python_function_legend || !s_python_function_ylim ||
                    !s_python_function_title || !s_python_function_axis ||
                    !s_python_function_xlabel || !s_python_function_ylabel ||
                    !s_python_function_xticks || !s_python_function_yticks ||
                    !s_python_function_xscale || !s_python_function_yscale ||
                    !s_python_function_grid || !s_python_function_xlim ||
                    !s_python_function_ion || !s_python_function_ginput ||
                    !s_python_function_save || !s_python_function_clf ||
                    !s_python_function_annotate || !s_python_function_errorbar ||
                    !s_python_function_errorbar || !s_python_function_tight_layout ||
                    !s_python_function_stem || !s_python_function_xkcd ||
                    !s_python_function_text || !s_python_function_suptitle ||
                    !s_python_function_bar || !s_python_function_subplots_adjust ||
                    !s_python_function_spy || !s_python_function_imshow) {
                    throw std::runtime_error("Couldn't find required function!");
                }

                if (!PyFunction_Check(s_python_function_show) ||
                    !PyFunction_Check(s_python_function_close) ||
                    !PyFunction_Check(s_python_function_draw) ||
                    !PyFunction_Check(s_python_function_pause) ||
                    !PyFunction_Check(s_python_function_figure) ||
                    !PyFunction_Check(s_python_function_fignum_exists) ||
                    !PyFunction_Check(s_python_function_plot) ||
                    !PyFunction_Check(s_python_function_quiver) ||
                    !PyFunction_Check(s_python_function_contour) ||
                    !PyFunction_Check(s_python_function_semilogx) ||
                    !PyFunction_Check(s_python_function_semilogy) ||
                    !PyFunction_Check(s_python_function_loglog) ||
                    !PyFunction_Check(s_python_function_fill) ||
                    !PyFunction_Check(s_python_function_fill_between) ||
                    !PyFunction_Check(s_python_function_spy) ||
                    !PyFunction_Check(s_python_function_subplot) ||
                    !PyFunction_Check(s_python_function_legend) ||
                    !PyFunction_Check(s_python_function_annotate) ||
                    !PyFunction_Check(s_python_function_ylim) ||
                    !PyFunction_Check(s_python_function_title) ||
                    !PyFunction_Check(s_python_function_axis) ||
                    !PyFunction_Check(s_python_function_xlabel) ||
                    !PyFunction_Check(s_python_function_ylabel) ||
                    !PyFunction_Check(s_python_function_xticks) ||
                    !PyFunction_Check(s_python_function_yticks) ||
                    !PyFunction_Check(s_python_function_xscale) ||
                    !PyFunction_Check(s_python_function_yscale) ||
                    !PyFunction_Check(s_python_function_grid) ||
                    !PyFunction_Check(s_python_function_xlim) ||
                    !PyFunction_Check(s_python_function_ion) ||
                    !PyFunction_Check(s_python_function_ginput) ||
                    !PyFunction_Check(s_python_function_save) ||
                    !PyFunction_Check(s_python_function_clf) ||
                    !PyFunction_Check(s_python_function_tight_layout) ||
                    !PyFunction_Check(s_python_function_errorbar) ||
                    !PyFunction_Check(s_python_function_stem) ||
                    !PyFunction_Check(s_python_function_xkcd) ||
                    !PyFunction_Check(s_python_function_text) ||
                    !PyFunction_Check(s_python_function_suptitle) ||
                    !PyFunction_Check(s_python_function_bar) ||
                    !PyFunction_Check(s_python_function_subplots_adjust) ||
                    !PyFunction_Check(s_python_function_imshow) ||
                    !PyFunction_Check(s_python_function_colorbar)
                    ) {
                    throw std::runtime_error(
                        "Python object is unexpectedly not a PyFunction.");
                }

                s_python_empty_tuple = PyTuple_New(0);
            }

            ~_interpreter() { Py_Finalize(); }
        };

    } // end namespace detail

    // must be called before the first regular call to matplotlib to have any effect
    inline void backend(const std::string& name) { detail::s_backend = name; }

    inline bool annotate(std::string annotation, double x, double y) {
        detail::_interpreter::get();

        PyObject* xy = PyTuple_New(2);
        PyObject* str = PyString_FromString(annotation.c_str());

        PyTuple_SetItem(xy, 0, PyFloat_FromDouble(x));
        PyTuple_SetItem(xy, 1, PyFloat_FromDouble(y));

        PyObject* kwargs = PyDict_New();
        PyDict_SetItemString(kwargs, "xy", xy);

        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, str);

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_annotate, args, kwargs);

        Py_DECREF(args);
        Py_DECREF(kwargs);

        if (res)
            Py_DECREF(res);

        return res;
    }

#ifndef WITHOUT_NUMPY
    // Type selector for numpy array conversion
    template <typename T> struct select_npy_type {
        const static NPY_TYPES type = NPY_NOTYPE;
    }; // Default
    template <> struct select_npy_type<double> {
        const static NPY_TYPES type = NPY_DOUBLE;
    };
    template <> struct select_npy_type<float> {
        const static NPY_TYPES type = NPY_FLOAT;
    };
    template <> struct select_npy_type<bool> {
        const static NPY_TYPES type = NPY_BOOL;
    };
    template <> struct select_npy_type<int8_t> {
        const static NPY_TYPES type = NPY_INT8;
    };
    template <> struct select_npy_type<int16_t> {
        const static NPY_TYPES type = NPY_SHORT;
    };
    template <> struct select_npy_type<int32_t> {
        const static NPY_TYPES type = NPY_INT;
    };
    template <> struct select_npy_type<int64_t> {
        const static NPY_TYPES type = NPY_INT64;
    };
    template <> struct select_npy_type<uint8_t> {
        const static NPY_TYPES type = NPY_UINT8;
    };
    template <> struct select_npy_type<uint16_t> {
        const static NPY_TYPES type = NPY_USHORT;
    };
    template <> struct select_npy_type<uint32_t> {
        const static NPY_TYPES type = NPY_ULONG;
    };
    template <> struct select_npy_type<uint64_t> {
        const static NPY_TYPES type = NPY_UINT64;
    };

    template <typename Vector> PyObject* get_array(const Vector& v) {
        detail::_interpreter::get(); // interpreter needs to be initialized for the
        // numpy commands to work
// both Eigen::Matrix<..> and std::vector<..> have the member value_type
        NPY_TYPES type = select_npy_type<typename Vector::value_type>::type;
        if (type == NPY_NOTYPE) {
            std::vector<double> vd(v.size());
            npy_intp vsize = v.size();
            // Eigen Vectors do not support begin/end() in the currently stable version
            // this can be changed once Eigen 3.4. is released
            // data() returns a pointer to the storage of the first element. If the
            // vector is modified afterwards, it may be rendedered invalid.
            // Note, that this is not an issue since get_array() is called by a
            // plot command using the instantaneous state of the vector.
            // The pointer is not reused. If the vector is plotted anew, data() is
            // called again and again get's the current, valid storage location.
            std::copy(v.data(), v.data() + v.size(), vd.begin());
            PyObject* varray =
                PyArray_SimpleNewFromData(1, &vsize, NPY_DOUBLE, (void*)(vd.data()));
            return varray;
        }

        npy_intp vsize = v.size();
        PyObject* varray =
            PyArray_SimpleNewFromData(1, &vsize, type, (void*)(v.data()));
        return varray;
    }

    // specialized get_2darray function for nested std::vectors
    template <typename Numeric>
    PyObject* get_2darray(const std::vector<::std::vector<Numeric>>& v) {
        detail::_interpreter::get(); // interpreter needs to be initialized for the
        // numpy commands to work
        if (v.size() < 1)
            throw std::runtime_error("get_2darray v too small");

        npy_intp vsize[2] = { static_cast<npy_intp>(v.size()),
                             static_cast<npy_intp>(v[0].size()) };

        PyArrayObject* varray =
            (PyArrayObject*)PyArray_SimpleNew(2, vsize, NPY_DOUBLE);

        double* vd_begin = static_cast<double*>(PyArray_DATA(varray));

        for (const ::std::vector<Numeric>& v_row : v) {
            if (v_row.size() != static_cast<size_t>(vsize[1]))
                throw std::runtime_error("mismatched array size");
            std::copy(v_row.begin(), v_row.end(), vd_begin);
            vd_begin += vsize[1];
        }

        return reinterpret_cast<PyObject*>(varray);
    }

    // suitable for more general matrices (especially Eigen matrices)
    template <typename Matrix> PyObject* get_2darray(const Matrix& A) {
        detail::_interpreter::get(); // interpreter needs to be initialized for the
        // numpy commands to work
        if (A.size() < 1)
            throw std::runtime_error("get_2darray A too small");

        npy_intp vsize[2] = { static_cast<npy_intp>(A.rows()),
                             static_cast<npy_intp>(A.cols()) };

        PyArrayObject* varray =
            (PyArrayObject*)PyArray_SimpleNew(2, vsize, NPY_DOUBLE);

        double* vd_begin = static_cast<double*>(PyArray_DATA(varray));

        for (std::size_t i = 0; i < A.rows(); ++i) {
            for (std::size_t j = 0; j < A.cols(); ++j) {
                *(vd_begin + i * A.cols() + j) = A(i, j);
            }
        }

        return reinterpret_cast<PyObject*>(varray);
    }

#else // fallback if we don't have numpy: copy every element of the given vector

    template <typename Vector> PyObject* get_array(const Vector& v) {
        detail::_interpreter::get();
        PyObject* list = PyList_New(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            PyList_SetItem(list, i, PyFloat_FromDouble(v.at(i)));
        }
        return list;
    }

#endif // WITHOUT_NUMPY

    namespace detail {
        // @brief Since most of the plot commands require the exact same usage apart
        //        from the call to the correct Python function, we encapsulate this
        // @param pyfunc The matplotlib function to be called with the given arguments
        // @param x The x vector
        // @param y The y vector
        // @param s The formatting string for colour, marker and linestyle
        // @param keywords Additional keywords, such as label
        // @return true if plot was successful, false otherwise
        template <typename VectorX, typename VectorY>
        bool plot_base(PyObject* const pyfunc, const VectorX& x, const VectorY& y,
            const std::string& s = "",
            const std::map<std::string, std::string>& keywords = {}) {
            assert(x.size() == y.size());

            PyObject* xarray = get_array(x);
            PyObject* yarray = get_array(y);

            PyObject* pystring = PyString_FromString(s.c_str());

            PyObject* plot_args = PyTuple_New(3);
            PyTuple_SetItem(plot_args, 0, xarray);
            PyTuple_SetItem(plot_args, 1, yarray);
            PyTuple_SetItem(plot_args, 2, pystring);

            PyObject* kwargs = PyDict_New();
            for (auto const& item : keywords) {
                PyDict_SetItemString(kwargs, item.first.c_str(),
                    PyString_FromString(item.second.c_str()));
            }

            PyObject* res = PyObject_Call(pyfunc, plot_args, kwargs);

            Py_DECREF(plot_args);
            Py_DECREF(kwargs);
            if (res)
                Py_DECREF(res);

            return res;
        }

    } // namespace detail

    // @brief standard plot function supporting the args (x, y, s, keywords)
    // @param x The x vector
    // @param y The y vector
    // @param s The formatting string
    // @param keywords Additional keywords
    // @return true, if successful, false otherwise
    template <typename VectorX, typename VectorY>
    bool plot(const VectorX& x, const VectorY& y, const std::string& s = "",
        const std::map<std::string, std::string>& keywords = {}) {
        return detail::plot_base(detail::_interpreter::get().s_python_function_plot,
            x, y, s, keywords);
    }

    // @brief standard plot function without formatting string, needed if
    //        keywords are given but formatting string is not
    template <typename VectorX, typename VectorY>
    bool plot(const VectorX& x, const VectorY& y,
        const std::map<std::string, std::string>& keywords) {
        return plot(x, y, "", keywords);
    }

    // @brief standard plot function if x data is not specified
    template <typename VectorY = std::vector<double>>
    bool plot(const VectorY& y, const std::string& format = "",
        const std::map<std::string, std::string>& keywords = {}) {
        // Note: cannot be an unsigned type for some reason, yields an overflow
        // problem..
        std::vector<int> x(y.size());
        for (int i = 0; i < x.size(); ++i)
            x.at(i) = i;

        return plot(x, y, format);
    }

    // @brief standard plot function if x data is not specified and the formatting
    //        string is missing
    template <typename VectorY = std::vector<double>>
    bool plot(const VectorY& y,
        const std::map<std::string, std::string>& keywords) {
        std::vector<int> x(y.size());
        for (int i = 0; i < x.size(); ++i)
            x.at(i) = i;

        return plot(x, y, "", keywords);
    }

    // @brief loglog plot function, see `plot` for more detail
    template <typename VectorX, typename VectorY>
    bool loglog(const VectorX& x, const VectorY& y, const std::string& s = "",
        const std::map<std::string, std::string>& keywords = {}) {
        return detail::plot_base(detail::_interpreter::get().s_python_function_loglog,
            x, y, s, keywords);
    }

    template <typename VectorX, typename VectorY>
    bool loglog(const VectorX& x, const VectorY& y,
        const std::map<std::string, std::string>& keywords) {
        return loglog(x, y, "", keywords);
    }

    template <typename VectorY = std::vector<double>>
    bool loglog(const VectorY& y, const std::string& s = "",
        const std::map<std::string, std::string>& keywords = {}) {
        std::vector<std::size_t> x(y.size());
        for (std::size_t i = 0; i < x.size(); ++i)
            x.at(i) = i;

        return loglog(x, y, s, keywords);
    }

    template <typename VectorY = std::vector<double>>
    bool loglog(const VectorY& y,
        const std::map<std::string, std::string>& keywords) {
        std::vector<std::size_t> x(y.size());
        for (std::size_t i = 0; i < x.size(); ++i)
            x.at(i) = i;

        return loglog(x, y, "", keywords);
    }

    // @brief semilogx plot function, see `plot` for more detail
    template <typename VectorX, typename VectorY>
    bool semilogx(const VectorX& x, const VectorY& y, const std::string& s = "",
        const std::map<std::string, std::string>& keywords = {}) {
        return detail::plot_base(
            detail::_interpreter::get().s_python_function_semilogx, x, y, s,
            keywords);
    }

    template <typename VectorX, typename VectorY>
    bool semilogx(const VectorX& x, const VectorY& y,
        const std::map<std::string, std::string>& keywords) {
        return semilogx(x, y, "", keywords);
    }

    template <typename VectorY = std::vector<double>>
    bool semilogx(const VectorY& y, const std::string& s = "",
        const std::map<std::string, std::string>& keywords = {}) {
        std::vector<std::size_t> x(y.size());
        for (std::size_t i = 0; i < x.size(); ++i)
            x.at(i) = i;

        return semilogx(x, y, s, keywords);
    }

    template <typename VectorY = std::vector<double>>
    bool semilogx(const VectorY& y,
        const std::map<std::string, std::string>& keywords) {
        std::vector<std::size_t> x(y.size());
        for (std::size_t i = 0; i < x.size(); ++i)
            x.at(i) = i;

        return semilogx(x, y, "", keywords);
    }

    // @brief semilogy plot function, see `plot` for more detail
    template <typename VectorX, typename VectorY>
    bool semilogy(const VectorX& x, const VectorY& y, const std::string& s = "",
        const std::map<std::string, std::string>& keywords = {}) {
        return detail::plot_base(
            detail::_interpreter::get().s_python_function_semilogy, x, y, s,
            keywords);
    }

    template <typename VectorX, typename VectorY>
    bool semilogy(const VectorX& x, const VectorY& y,
        const std::map<std::string, std::string>& keywords) {
        return semilogy(x, y, "", keywords);
    }

    template <typename VectorY = std::vector<double>>
    bool semilogy(const VectorY& y, const std::string& s = "",
        const std::map<std::string, std::string>& keywords = {}) {
        std::vector<std::size_t> x(y.size());
        for (std::size_t i = 0; i < x.size(); ++i)
            x.at(i) = i;

        return semilogy(x, y, s, keywords);
    }

    template <typename VectorY = std::vector<double>>
    bool semilogy(const VectorY& y,
        const std::map<std::string, std::string>& keywords) {
        std::vector<std::size_t> x(y.size());
        for (std::size_t i = 0; i < x.size(); ++i)
            x.at(i) = i;

        return semilogy(x, y, "", keywords);
    }

    template <typename Matrix>
    void imshow(const Matrix& X, const std::map<std::string, std::string>& keywords = {}) {
        PyObject* Xarray = get_2darray(X);

        PyObject* kwargs = PyDict_New();
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* plot_args = PyTuple_New(1);
        PyTuple_SetItem(plot_args, 0, Xarray);

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_imshow, plot_args, kwargs);

        Py_DECREF(plot_args);
        Py_DECREF(kwargs);
        if (res)
            Py_DECREF(res);
    }

    // @brief Add the colorbar
    void colorbar() {
        PyObject* res =
            PyObject_CallObject(detail::_interpreter::get().s_python_function_colorbar,
                detail::_interpreter::get().s_python_empty_tuple);
        if (!res)
            throw std::runtime_error("Call to colorbar() failed.");
    }

    // @brief plot_surface for datapoints (x_ij, y_ij, z_ij) with i,j = 0..n
    // @param x The x values of the datapoints in a matrix
    // @param y The y values of the datapoints in a matrix
    // @param z The function value of the datapoints in a matrix
    // @param keywords Additional keywords
    template <typename Matrix>
    void plot_surface(const Matrix& x, const Matrix& y, const Matrix& z,
        const std::map<std::string, std::string>& keywords =
        std::map<std::string, std::string>()) {
        // We lazily load the modules here the first time this function is called
        // because I'm not sure that we can assume "matplotlib installed" implies
        // "mpl_toolkits installed" on all platforms, and we don't want to require
        // it for people who don't need 3d plots.
        static PyObject* mpl_toolkitsmod = nullptr, * axis3dmod = nullptr;
        if (!mpl_toolkitsmod) {
            detail::_interpreter::get();

            PyObject* mpl_toolkits = PyString_FromString("mpl_toolkits");
            PyObject* axis3d = PyString_FromString("mpl_toolkits.mplot3d");
            if (!mpl_toolkits || !axis3d) {
                throw std::runtime_error("couldnt create string");
            }

            mpl_toolkitsmod = PyImport_Import(mpl_toolkits);
            Py_DECREF(mpl_toolkits);
            if (!mpl_toolkitsmod) {
                throw std::runtime_error("Error loading module mpl_toolkits!");
            }

            axis3dmod = PyImport_Import(axis3d);
            Py_DECREF(axis3d);
            if (!axis3dmod) {
                throw std::runtime_error("Error loading module mpl_toolkits.mplot3d!");
            }
        }

        assert(x.size() == y.size());
        assert(y.size() == z.size());

        // using numpy arrays
        PyObject* xarray = get_2darray(x);
        PyObject* yarray = get_2darray(y);
        PyObject* zarray = get_2darray(z);

        // construct positional args
        PyObject* args = PyTuple_New(3);
        PyTuple_SetItem(args, 0, xarray);
        PyTuple_SetItem(args, 1, yarray);
        PyTuple_SetItem(args, 2, zarray);

        // Build up the kw args.
        PyObject* kwargs = PyDict_New();
        PyDict_SetItemString(kwargs, "rstride", PyInt_FromLong(1));
        PyDict_SetItemString(kwargs, "cstride", PyInt_FromLong(1));

        PyObject* python_colormap_coolwarm = PyObject_GetAttrString(
            detail::_interpreter::get().s_python_colormap, "coolwarm");

        PyDict_SetItemString(kwargs, "cmap", python_colormap_coolwarm);

        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyString_FromString(it->second.c_str()));
        }

        PyObject* fig =
            PyObject_CallObject(detail::_interpreter::get().s_python_function_figure,
                detail::_interpreter::get().s_python_empty_tuple);
        if (!fig)
            throw std::runtime_error("Call to figure() failed.");

        PyObject* gca_kwargs = PyDict_New();
        PyDict_SetItemString(gca_kwargs, "projection", PyString_FromString("3d"));

        PyObject* gca = PyObject_GetAttrString(fig, "gca");
        if (!gca)
            throw std::runtime_error("No gca");
        Py_INCREF(gca);
        PyObject* axis = PyObject_Call(
            gca, detail::_interpreter::get().s_python_empty_tuple, gca_kwargs);

        if (!axis)
            throw std::runtime_error("No axis");
        Py_INCREF(axis);

        Py_DECREF(gca);
        Py_DECREF(gca_kwargs);

        PyObject* plot_surface = PyObject_GetAttrString(axis, "plot_surface");
        if (!plot_surface)
            throw std::runtime_error("No surface");
        Py_INCREF(plot_surface);
        PyObject* res = PyObject_Call(plot_surface, args, kwargs);
        if (!res)
            throw std::runtime_error("failed surface");
        Py_DECREF(plot_surface);

        Py_DECREF(axis);
        Py_DECREF(args);
        Py_DECREF(kwargs);
        if (res)
            Py_DECREF(res);
    }

    // @brief plot_surface for datapoints (x_ij, y_ij, z_ij) with i,j = 0..n
    // @param x The x values of the datapoints in a matrix
    // @param y The y values of the datapoints in a matrix
    // @param z The function value of the datapoints in a matrix
    // @param keywords Additional keywords
    template <typename Matrix>
    void contour(const Matrix& x, const Matrix& y, const Matrix& z,
        const std::map<std::string, std::string>& keywords = {}) {
        detail::_interpreter::get();

        // using numpy arrays
        PyObject* xarray = get_2darray(x);
        PyObject* yarray = get_2darray(y);
        PyObject* zarray = get_2darray(z);

        // construct positional args
        PyObject* args = PyTuple_New(3);
        PyTuple_SetItem(args, 0, xarray);
        PyTuple_SetItem(args, 1, yarray);
        PyTuple_SetItem(args, 2, zarray);

        // Build up the kw args.
        PyObject* kwargs = PyDict_New();

        PyObject* python_colormap_coolwarm = PyObject_GetAttrString(
            detail::_interpreter::get().s_python_colormap, "coolwarm");

        PyDict_SetItemString(kwargs, "cmap", python_colormap_coolwarm);

        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyString_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_contour, args, kwargs);
        if (!res)
            throw std::runtime_error("failed surface");

        Py_DECREF(args);
        Py_DECREF(kwargs);
        if (res)
            Py_DECREF(res);
    }

    template <typename Numeric>
    bool stem(const std::vector<Numeric>& x, const std::vector<Numeric>& y,
        const std::map<std::string, std::string>& keywords) {
        assert(x.size() == y.size());

        // using numpy arrays
        PyObject* xarray = get_array(x);
        PyObject* yarray = get_array(y);

        // construct positional args
        PyObject* args = PyTuple_New(2);
        PyTuple_SetItem(args, 0, xarray);
        PyTuple_SetItem(args, 1, yarray);

        // construct keyword args
        PyObject* kwargs = PyDict_New();
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyString_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_stem, args, kwargs);

        Py_DECREF(args);
        Py_DECREF(kwargs);
        if (res)
            Py_DECREF(res);

        return res;
    }

    template <typename Numeric>
    bool fill(const std::vector<Numeric>& x, const std::vector<Numeric>& y,
        const std::map<std::string, std::string>& keywords) {
        assert(x.size() == y.size());

        // using numpy arrays
        PyObject* xarray = get_array(x);
        PyObject* yarray = get_array(y);

        // construct positional args
        PyObject* args = PyTuple_New(2);
        PyTuple_SetItem(args, 0, xarray);
        PyTuple_SetItem(args, 1, yarray);

        // construct keyword args
        PyObject* kwargs = PyDict_New();
        for (auto it = keywords.begin(); it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_fill, args, kwargs);

        Py_DECREF(args);
        Py_DECREF(kwargs);

        if (res)
            Py_DECREF(res);

        return res;
    }

    template <typename Numeric>
    bool fill_between(const std::vector<Numeric>& x, const std::vector<Numeric>& y1,
        const std::vector<Numeric>& y2,
        const std::map<std::string, std::string>& keywords) {
        assert(x.size() == y1.size());
        assert(x.size() == y2.size());

        // using numpy arrays
        PyObject* xarray = get_array(x);
        PyObject* y1array = get_array(y1);
        PyObject* y2array = get_array(y2);

        // construct positional args
        PyObject* args = PyTuple_New(3);
        PyTuple_SetItem(args, 0, xarray);
        PyTuple_SetItem(args, 1, y1array);
        PyTuple_SetItem(args, 2, y2array);

        // construct keyword args
        PyObject* kwargs = PyDict_New();
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_fill_between, args, kwargs);

        Py_DECREF(args);
        Py_DECREF(kwargs);
        if (res)
            Py_DECREF(res);

        return res;
    }

    template <typename VectorY>
    bool hist(const VectorY& y, long bins = 10, std::string color = "b",
        double alpha = 1.0, bool cumulative = false) {

        PyObject* yarray = get_array(y);

        PyObject* kwargs = PyDict_New();
        PyDict_SetItemString(kwargs, "bins", PyLong_FromLong(bins));
        PyDict_SetItemString(kwargs, "color", PyString_FromString(color.c_str()));
        PyDict_SetItemString(kwargs, "alpha", PyFloat_FromDouble(alpha));
        PyDict_SetItemString(kwargs, "cumulative", cumulative ? Py_True : Py_False);

        PyObject* plot_args = PyTuple_New(1);

        PyTuple_SetItem(plot_args, 0, yarray);

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_hist, plot_args, kwargs);

        Py_DECREF(plot_args);
        Py_DECREF(kwargs);
        if (res)
            Py_DECREF(res);

        return res;
    }

    // @brief Scatter plot
    // @param x x-coordinates of the 2d points
    // @param y y-coordinates of the 2d points
    // @param s the marker size in points**2
    // @param keywords Additional keywords
    template <typename VectorX, typename VectorY>
    bool scatter(const VectorX& x, const VectorY& y, const double s = 1.0,
        const std::map<std::string, std::string> keywords = {}) {
        assert(x.size() == y.size());

        PyObject* xarray = get_array(x);
        PyObject* yarray = get_array(y);

        PyObject* kwargs = PyDict_New();
        PyDict_SetItemString(kwargs, "s", PyFloat_FromDouble(s));
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* plot_args = PyTuple_New(2);
        PyTuple_SetItem(plot_args, 0, xarray);
        PyTuple_SetItem(plot_args, 1, yarray);

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_scatter, plot_args, kwargs);

        Py_DECREF(plot_args);
        Py_DECREF(kwargs);
        if (res)
            Py_DECREF(res);

        return res;
    }

    template <typename VectorX, typename VectorY>
    bool scatter(const VectorX& x, const VectorY& y,
        const std::map<std::string, std::string>& keywords) {
        return scatter(x, y, 1.0, keywords);
    }

    // @brief Spy plot
    // @param A the matrix
    // @param precision Plot all elements above `|precision|`
    // @param keywords Additional keywords
    template <typename Matrix>
    bool spy(const Matrix& A,
        const std::map<std::string, std::string>& keywords = {}) {
        PyObject* Aarray = get_2darray(A);

        PyObject* kwargs = PyDict_New();
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* plot_args = PyTuple_New(1);
        PyTuple_SetItem(plot_args, 0, Aarray);

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_spy, plot_args, kwargs);

        Py_DECREF(plot_args);
        Py_DECREF(kwargs);
        if (res)
            Py_DECREF(res);

        return res;
    }

    template <typename Numeric>
    bool bar(const std::vector<Numeric>& y, std::string ec = "black",
        std::string ls = "-", double lw = 1.0,
        const std::map<std::string, std::string>& keywords = {}) {
        PyObject* yarray = get_array(y);

        std::vector<int> x;
        for (int i = 0; i < y.size(); i++)
            x.push_back(i);

        PyObject* xarray = get_array(x);

        PyObject* kwargs = PyDict_New();

        PyDict_SetItemString(kwargs, "ec", PyString_FromString(ec.c_str()));
        PyDict_SetItemString(kwargs, "ls", PyString_FromString(ls.c_str()));
        PyDict_SetItemString(kwargs, "lw", PyFloat_FromDouble(lw));

        PyObject* plot_args = PyTuple_New(2);
        PyTuple_SetItem(plot_args, 0, xarray);
        PyTuple_SetItem(plot_args, 1, yarray);

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_bar, plot_args, kwargs);

        Py_DECREF(plot_args);
        Py_DECREF(kwargs);
        if (res)
            Py_DECREF(res);

        return res;
    }

    inline bool
        subplots_adjust(const std::map<std::string, double>& keywords = {}) {

        PyObject* kwargs = PyDict_New();
        for (std::map<std::string, double>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyFloat_FromDouble(it->second));
        }

        PyObject* plot_args = PyTuple_New(0);

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_subplots_adjust, plot_args,
            kwargs);

        Py_DECREF(plot_args);
        Py_DECREF(kwargs);
        if (res)
            Py_DECREF(res);

        return res;
    }

    template <typename Numeric>
    bool named_hist(std::string label, const std::vector<Numeric>& y,
        long bins = 10, std::string color = "b", double alpha = 1.0) {
        PyObject* yarray = get_array(y);

        PyObject* kwargs = PyDict_New();
        PyDict_SetItemString(kwargs, "label", PyString_FromString(label.c_str()));
        PyDict_SetItemString(kwargs, "bins", PyLong_FromLong(bins));
        PyDict_SetItemString(kwargs, "color", PyString_FromString(color.c_str()));
        PyDict_SetItemString(kwargs, "alpha", PyFloat_FromDouble(alpha));

        PyObject* plot_args = PyTuple_New(1);
        PyTuple_SetItem(plot_args, 0, yarray);

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_hist, plot_args, kwargs);

        Py_DECREF(plot_args);
        Py_DECREF(kwargs);
        if (res)
            Py_DECREF(res);

        return res;
    }

    template <typename NumericX, typename NumericY, typename NumericU,
        typename NumericW>
    bool quiver(const std::vector<NumericX>& x, const std::vector<NumericY>& y,
        const std::vector<NumericU>& u, const std::vector<NumericW>& w,
        const std::map<std::string, std::string>& keywords = {}) {
        assert(x.size() == y.size() && x.size() == u.size() && u.size() == w.size());

        PyObject* xarray = get_array(x);
        PyObject* yarray = get_array(y);
        PyObject* uarray = get_array(u);
        PyObject* warray = get_array(w);

        PyObject* plot_args = PyTuple_New(4);
        PyTuple_SetItem(plot_args, 0, xarray);
        PyTuple_SetItem(plot_args, 1, yarray);
        PyTuple_SetItem(plot_args, 2, uarray);
        PyTuple_SetItem(plot_args, 3, warray);

        // construct keyword args
        PyObject* kwargs = PyDict_New();
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_quiver, plot_args, kwargs);

        Py_DECREF(kwargs);
        Py_DECREF(plot_args);
        if (res)
            Py_DECREF(res);

        return res;
    }

    template <typename NumericY>
    void axhline(const NumericY y,
        const std::map<std::string, std::string> keywords = {}) {
        detail::_interpreter::get();

        PyObject* kwargs = PyDict_New();

        // add location
        PyDict_SetItemString(kwargs, "y", PyFloat_FromDouble(y));

        // add other keywords
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* res =
            PyObject_Call(detail::_interpreter::get().s_python_function_axhline,
                detail::_interpreter::get().s_python_empty_tuple, kwargs);

        Py_DECREF(kwargs);

        if (!res)
            throw std::runtime_error("Call to axhline() failed.");

        Py_DECREF(res);
    }

    template <typename NumericX>
    void axvline(const NumericX x,
        const std::map<std::string, std::string> keywords = {}) {
        detail::_interpreter::get();

        PyObject* kwargs = PyDict_New();

        // add location
        PyDict_SetItemString(kwargs, "x", PyFloat_FromDouble(x));

        // add other keywords
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* res =
            PyObject_Call(detail::_interpreter::get().s_python_function_axvline,
                detail::_interpreter::get().s_python_empty_tuple, kwargs);

        Py_DECREF(kwargs);

        if (!res)
            throw std::runtime_error("Call to axvline() failed.");

        Py_DECREF(res);
    }

    template <typename NumericX, typename NumericY>
    bool stem(const std::vector<NumericX>& x, const std::vector<NumericY>& y,
        const std::string& s = "") {
        assert(x.size() == y.size());

        PyObject* xarray = get_array(x);
        PyObject* yarray = get_array(y);

        PyObject* pystring = PyString_FromString(s.c_str());

        PyObject* plot_args = PyTuple_New(3);
        PyTuple_SetItem(plot_args, 0, xarray);
        PyTuple_SetItem(plot_args, 1, yarray);
        PyTuple_SetItem(plot_args, 2, pystring);

        PyObject* res = PyObject_CallObject(
            detail::_interpreter::get().s_python_function_stem, plot_args);

        Py_DECREF(plot_args);
        if (res)
            Py_DECREF(res);

        return res;
    }

    template <typename VectorX, typename VectorY>
    bool errorbar(const VectorX& x, const VectorY& y, const VectorY& yerr,
        const std::map<std::string, std::string>& keywords = {}) {
        assert(x.size() == y.size());

        PyObject* xarray = get_array(x);
        PyObject* yarray = get_array(y);
        PyObject* yerrarray = get_array(yerr);

        // construct keyword args
        PyObject* kwargs = PyDict_New();
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyString_FromString(it->second.c_str()));
        }

        PyDict_SetItemString(kwargs, "yerr", yerrarray);

        PyObject* plot_args = PyTuple_New(2);
        PyTuple_SetItem(plot_args, 0, xarray);
        PyTuple_SetItem(plot_args, 1, yarray);

        PyObject* res =
            PyObject_Call(detail::_interpreter::get().s_python_function_errorbar,
                plot_args, kwargs);

        Py_DECREF(kwargs);
        Py_DECREF(plot_args);

        if (res)
            Py_DECREF(res);
        else
            throw std::runtime_error("Call to errorbar() failed.");

        return res;
    }

    template <typename Numeric>
    bool stem(const std::vector<Numeric>& y, const std::string& format = "") {
        std::vector<Numeric> x(y.size());
        for (size_t i = 0; i < x.size(); ++i)
            x.at(i) = i;
        return stem(x, y, format);
    }

    template <typename Numeric>
    void text(Numeric x, Numeric y, const std::string& s = "") {
        detail::_interpreter::get();

        PyObject* args = PyTuple_New(3);
        PyTuple_SetItem(args, 0, PyFloat_FromDouble(x));
        PyTuple_SetItem(args, 1, PyFloat_FromDouble(y));
        PyTuple_SetItem(args, 2, PyString_FromString(s.c_str()));

        PyObject* res = PyObject_CallObject(
            detail::_interpreter::get().s_python_function_text, args);
        if (!res)
            throw std::runtime_error("Call to text() failed.");

        Py_DECREF(args);
        Py_DECREF(res);
    }

    inline long figure(long number = -1) {
        PyObject* res;
        if (number == -1)
            res = PyObject_CallObject(
                detail::_interpreter::get().s_python_function_figure,
                detail::_interpreter::get().s_python_empty_tuple);
        else {
            assert(number > 0);

            // Make sure interpreter is initialised
            detail::_interpreter::get();

            PyObject* args = PyTuple_New(1);
            PyTuple_SetItem(args, 0, PyLong_FromLong(number));
            res = PyObject_CallObject(
                detail::_interpreter::get().s_python_function_figure, args);
            Py_DECREF(args);
        }

        if (!res)
            throw std::runtime_error("Call to figure() failed.");

        PyObject* num = PyObject_GetAttrString(res, "number");
        if (!num)
            throw std::runtime_error("Could not get number attribute of figure object");
        const long figureNumber = PyLong_AsLong(num);

        Py_DECREF(num);
        Py_DECREF(res);

        return figureNumber;
    }

    inline bool fignum_exists(long number) {
        // Make sure interpreter is initialised
        detail::_interpreter::get();

        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, PyLong_FromLong(number));
        PyObject* res = PyObject_CallObject(
            detail::_interpreter::get().s_python_function_fignum_exists, args);
        if (!res)
            throw std::runtime_error("Call to fignum_exists() failed.");

        bool ret = PyObject_IsTrue(res);
        Py_DECREF(res);
        Py_DECREF(args);

        return ret;
    }

    inline void figure_size(size_t w, size_t h) {
        // Make sure interpreter is initialised
        detail::_interpreter::get();

        const size_t dpi = 100;
        PyObject* size = PyTuple_New(2);
        PyTuple_SetItem(size, 0, PyFloat_FromDouble((double)w / dpi));
        PyTuple_SetItem(size, 1, PyFloat_FromDouble((double)h / dpi));

        PyObject* kwargs = PyDict_New();
        PyDict_SetItemString(kwargs, "figsize", size);
        PyDict_SetItemString(kwargs, "dpi", PyLong_FromSize_t(dpi));

        PyObject* res =
            PyObject_Call(detail::_interpreter::get().s_python_function_figure,
                detail::_interpreter::get().s_python_empty_tuple, kwargs);

        Py_DECREF(kwargs);

        if (!res)
            throw std::runtime_error("Call to figure_size() failed.");
        Py_DECREF(res);
    }

    template <typename Vector = std::vector<double>>
    inline void legend(const std::string& loc = "best",
        const Vector& bbox_to_anchor = Vector(),
        const std::map<std::string, std::string>& keywords = {}) {
        detail::_interpreter::get();

        PyObject* kwargs = PyDict_New();

        // add location
        if (loc != "")
            PyDict_SetItemString(kwargs, "loc", PyString_FromString(loc.c_str()));

        // add bbox to anchor
        if (bbox_to_anchor.size() == 2 || bbox_to_anchor.size() == 4) {
            PyObject* bbox = get_array(bbox_to_anchor);
            PyDict_SetItemString(kwargs, "bbox_to_anchor", bbox);
        }

        // add other keywords
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyString_FromString(it->second.c_str()));
        }

        PyObject* res =
            PyObject_Call(detail::_interpreter::get().s_python_function_legend,
                detail::_interpreter::get().s_python_empty_tuple, kwargs);

        Py_DECREF(kwargs);

        if (!res)
            throw std::runtime_error("Call to legend() failed.");

        Py_DECREF(res);
    }

    template <typename Vector>
    inline void legend(const Vector& bbox_to_anchor,
        const std::map<std::string, std::string>& keywords = {}) {
        legend("", bbox_to_anchor, keywords);
    }

    inline void legend(const std::string& loc,
        const std::map<std::string, std::string>& keywords = {}) {
        legend(loc, std::vector<double>(), keywords);
    }

    // to support C-style strings we also need const char[], std::string only
    // does not capture calls of style legend("lower left")
    inline void legend(const char loc[],
        const std::map<std::string, std::string>& keywords = {}) {
        legend(loc, std::vector<double>(), keywords);
    }

    /*
    inline void legend(const std::string& loc,
                       const std::map<std::string, std::string>& keywords = {}) {
      legend(loc, std::vector<double>(), keywords);
    }

    inline void legend(const std::map<std::string, std::string>& keywords) {
      legend("", std::vector<double>(), keywords);
    }
    */

    template <typename Numeric> void ylim(const Numeric bottom, const Numeric top) {
        detail::_interpreter::get();

        PyObject* list = PyList_New(2);
        PyList_SetItem(list, 0, PyFloat_FromDouble(bottom));
        PyList_SetItem(list, 1, PyFloat_FromDouble(top));

        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, list);

        PyObject* res = PyObject_CallObject(
            detail::_interpreter::get().s_python_function_ylim, args);
        if (!res)
            throw std::runtime_error("Call to ylim() failed.");

        Py_DECREF(args);
        Py_DECREF(res);
    }

    template <typename Numeric> void xlim(const Numeric left, const Numeric right) {
        detail::_interpreter::get();

        PyObject* list = PyList_New(2);
        PyList_SetItem(list, 0, PyFloat_FromDouble(left));
        PyList_SetItem(list, 1, PyFloat_FromDouble(right));

        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, list);

        PyObject* res = PyObject_CallObject(
            detail::_interpreter::get().s_python_function_xlim, args);
        if (!res)
            throw std::runtime_error("Call to xlim() failed.");

        Py_DECREF(args);
        Py_DECREF(res);
    }

    inline double* xlim() {
        detail::_interpreter::get();

        PyObject* args = PyTuple_New(0);
        PyObject* res = PyObject_CallObject(
            detail::_interpreter::get().s_python_function_xlim, args);
        PyObject* left = PyTuple_GetItem(res, 0);
        PyObject* right = PyTuple_GetItem(res, 1);

        double* arr = new double[2];
        arr[0] = PyFloat_AsDouble(left);
        arr[1] = PyFloat_AsDouble(right);

        if (!res)
            throw std::runtime_error("Call to xlim() failed.");

        Py_DECREF(res);
        return arr;
    }

    inline double* ylim() {
        detail::_interpreter::get();

        PyObject* args = PyTuple_New(0);
        PyObject* res = PyObject_CallObject(
            detail::_interpreter::get().s_python_function_ylim, args);
        PyObject* bottom = PyTuple_GetItem(res, 0);
        PyObject* top = PyTuple_GetItem(res, 1);

        double* arr = new double[2];
        arr[0] = PyFloat_AsDouble(bottom);
        arr[1] = PyFloat_AsDouble(top);

        if (!res)
            throw std::runtime_error("Call to ylim() failed.");

        Py_DECREF(res);
        return arr;
    }

    template <typename Numeric>
    inline void xticks(const std::vector<Numeric>& ticks,
        const std::vector<std::string>& labels = {},
        const std::map<std::string, std::string>& keywords = {}) {
        assert(labels.size() == 0 || ticks.size() == labels.size());

        // using numpy array
        PyObject* ticksarray = get_array(ticks);

        PyObject* args;
        if (labels.size() == 0) {
            // construct positional args
            args = PyTuple_New(1);
            PyTuple_SetItem(args, 0, ticksarray);
        }
        else {
            // make tuple of tick labels
            PyObject* labelstuple = PyTuple_New(labels.size());
            for (size_t i = 0; i < labels.size(); i++)
                PyTuple_SetItem(labelstuple, i, PyUnicode_FromString(labels[i].c_str()));

            // construct positional args
            args = PyTuple_New(2);
            PyTuple_SetItem(args, 0, ticksarray);
            PyTuple_SetItem(args, 1, labelstuple);
        }

        // construct keyword args
        PyObject* kwargs = PyDict_New();
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyString_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_xticks, args, kwargs);

        Py_DECREF(args);
        Py_DECREF(kwargs);
        if (!res)
            throw std::runtime_error("Call to xticks() failed");

        Py_DECREF(res);
    }

    template <typename Numeric>
    inline void xticks(const std::vector<Numeric>& ticks,
        const std::map<std::string, std::string>& keywords) {
        xticks(ticks, {}, keywords);
    }

    template <typename Numeric>
    inline void yticks(const std::vector<Numeric>& ticks,
        const std::vector<std::string>& labels = {},
        const std::map<std::string, std::string>& keywords = {}) {
        assert(labels.size() == 0 || ticks.size() == labels.size());

        // using numpy array
        PyObject* ticksarray = get_array(ticks);

        PyObject* args;
        if (labels.size() == 0) {
            // construct positional args
            args = PyTuple_New(1);
            PyTuple_SetItem(args, 0, ticksarray);
        }
        else {
            // make tuple of tick labels
            PyObject* labelstuple = PyTuple_New(labels.size());
            for (size_t i = 0; i < labels.size(); i++)
                PyTuple_SetItem(labelstuple, i, PyUnicode_FromString(labels[i].c_str()));

            // construct positional args
            args = PyTuple_New(2);
            PyTuple_SetItem(args, 0, ticksarray);
            PyTuple_SetItem(args, 1, labelstuple);
        }

        // construct keyword args
        PyObject* kwargs = PyDict_New();
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyString_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_yticks, args, kwargs);

        Py_DECREF(args);
        Py_DECREF(kwargs);
        if (!res)
            throw std::runtime_error("Call to yticks() failed");

        Py_DECREF(res);
    }

    template <typename Numeric>
    inline void yticks(const std::vector<Numeric>& ticks,
        const std::map<std::string, std::string>& keywords) {
        yticks(ticks, {}, keywords);
    }

    inline void subplot(long nrows, long ncols, long plot_number) {
        detail::_interpreter::get();

        // construct positional args
        PyObject* args = PyTuple_New(3);
        PyTuple_SetItem(args, 0, PyFloat_FromDouble(nrows));
        PyTuple_SetItem(args, 1, PyFloat_FromDouble(ncols));
        PyTuple_SetItem(args, 2, PyFloat_FromDouble(plot_number));

        PyObject* res = PyObject_CallObject(
            detail::_interpreter::get().s_python_function_subplot, args);
        if (!res)
            throw std::runtime_error("Call to subplot() failed.");

        Py_DECREF(args);
        Py_DECREF(res);
    }

    inline void title(const std::string& titlestr,
        const std::map<std::string, std::string>& keywords = {}) {
        detail::_interpreter::get();

        PyObject* pytitlestr = PyString_FromString(titlestr.c_str());
        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, pytitlestr);

        PyObject* kwargs = PyDict_New();
        for (auto it = keywords.begin(); it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_title, args, kwargs);
        if (!res)
            throw std::runtime_error("Call to title() failed.");

        Py_DECREF(args);
        Py_DECREF(kwargs);
        Py_DECREF(res);
    }

    inline void suptitle(const std::string& suptitlestr,
        const std::map<std::string, std::string>& keywords = {}) {
        detail::_interpreter::get();

        PyObject* pysuptitlestr = PyString_FromString(suptitlestr.c_str());
        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, pysuptitlestr);

        PyObject* kwargs = PyDict_New();
        for (auto it = keywords.begin(); it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_suptitle, args, kwargs);
        if (!res)
            throw std::runtime_error("Call to suptitle() failed.");

        Py_DECREF(args);
        Py_DECREF(kwargs);
        Py_DECREF(res);
    }

    inline void axis(const std::string& option) {
        PyObject* str = PyString_FromString(option.c_str());
        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, str);

        PyObject* res = PyObject_CallObject(
            detail::_interpreter::get().s_python_function_axis, args);
        if (!res)
            throw std::runtime_error("Call to axis() failed.");

        Py_DECREF(args);
        Py_DECREF(res);
    }

    inline void xlabel(const std::string& str,
        const std::map<std::string, std::string>& keywords = {}) {
        PyObject* pystr = PyString_FromString(str.c_str());
        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, pystr);

        PyObject* kwargs = PyDict_New();
        for (auto it = keywords.begin(); it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_xlabel, args, kwargs);
        if (!res)
            throw std::runtime_error("Call to xlabel() failed.");

        Py_DECREF(args);
        Py_DECREF(kwargs);
        Py_DECREF(res);
    }

    inline void ylabel(const std::string& str,
        const std::map<std::string, std::string>& keywords = {}) {
        PyObject* pystr = PyString_FromString(str.c_str());
        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, pystr);

        PyObject* kwargs = PyDict_New();
        for (auto it = keywords.begin(); it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_ylabel, args, kwargs);
        if (!res)
            throw std::runtime_error("Call to ylabel() failed.");

        Py_DECREF(args);
        Py_DECREF(kwargs);
        Py_DECREF(res);
    }

    inline void grid(bool flag = true) {
        PyObject* pyflag = flag ? Py_True : Py_False;
        Py_INCREF(pyflag);

        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, pyflag);

        PyObject* res = PyObject_CallObject(
            detail::_interpreter::get().s_python_function_grid, args);
        if (!res)
            throw std::runtime_error("Call to grid() failed.");

        Py_DECREF(args);
        Py_DECREF(res);
    }

    inline void show(const bool block = true) {
        PyObject* res;
        if (block) {
            res =
                PyObject_CallObject(detail::_interpreter::get().s_python_function_show,
                    detail::_interpreter::get().s_python_empty_tuple);
        }
        else {
            PyObject* kwargs = PyDict_New();
            PyDict_SetItemString(kwargs, "block", Py_False);
            res =
                PyObject_Call(detail::_interpreter::get().s_python_function_show,
                    detail::_interpreter::get().s_python_empty_tuple, kwargs);
            Py_DECREF(kwargs);
        }

        if (!res)
            throw std::runtime_error("Call to show() failed.");

        Py_DECREF(res);
    }

    inline void close() {
        PyObject* res =
            PyObject_CallObject(detail::_interpreter::get().s_python_function_close,
                detail::_interpreter::get().s_python_empty_tuple);

        if (!res)
            throw std::runtime_error("Call to close() failed.");

        Py_DECREF(res);
    }

    inline void xkcd() {
        PyObject* res;
        PyObject* kwargs = PyDict_New();

        res = PyObject_Call(detail::_interpreter::get().s_python_function_xkcd,
            detail::_interpreter::get().s_python_empty_tuple, kwargs);

        Py_DECREF(kwargs);

        if (!res)
            throw std::runtime_error("Call to xkcd() failed.");

        Py_DECREF(res);
    }

    inline void draw() {
        PyObject* res =
            PyObject_CallObject(detail::_interpreter::get().s_python_function_draw,
                detail::_interpreter::get().s_python_empty_tuple);

        if (!res)
            throw std::runtime_error("Call to draw() failed.");

        Py_DECREF(res);
    }

    template <typename Numeric> inline void pause(Numeric interval) {
        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, PyFloat_FromDouble(interval));

        PyObject* res = PyObject_CallObject(
            detail::_interpreter::get().s_python_function_pause, args);
        if (!res)
            throw std::runtime_error("Call to pause() failed.");

        Py_DECREF(args);
        Py_DECREF(res);
    }

    inline void savefig(const std::string& filename,
        const std::map<std::string, std::string>& keywords = {}) {
        PyObject* pyfilename = PyString_FromString(filename.c_str());

        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, pyfilename);

        PyObject* kwargs = PyDict_New();
        for (auto it = keywords.begin(); it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_save, args, kwargs);
        if (!res)
            throw std::runtime_error("Call to save() failed.");

        Py_DECREF(kwargs);
        Py_DECREF(args);
        Py_DECREF(res);
    }

    inline void save(const std::string& filename) {
        std::cerr << "matplotlibcpp::save is deprecated, use savefig instead\n";
        matplotlibcpp::savefig(filename);
    }

    inline void clf() {
        PyObject* res =
            PyObject_CallObject(detail::_interpreter::get().s_python_function_clf,
                detail::_interpreter::get().s_python_empty_tuple);

        if (!res)
            throw std::runtime_error("Call to clf() failed.");

        Py_DECREF(res);
    }

    inline void ion() {
        PyObject* res =
            PyObject_CallObject(detail::_interpreter::get().s_python_function_ion,
                detail::_interpreter::get().s_python_empty_tuple);

        if (!res)
            throw std::runtime_error("Call to ion() failed.");

        Py_DECREF(res);
    }

    inline std::vector<std::array<double, 2>>
        ginput(const int numClicks = 1,
            const std::map<std::string, std::string>& keywords = {}) {
        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, PyLong_FromLong(numClicks));

        // construct keyword args
        PyObject* kwargs = PyDict_New();
        for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
            it != keywords.end(); ++it) {
            PyDict_SetItemString(kwargs, it->first.c_str(),
                PyUnicode_FromString(it->second.c_str()));
        }

        PyObject* res = PyObject_Call(
            detail::_interpreter::get().s_python_function_ginput, args, kwargs);

        Py_DECREF(kwargs);
        Py_DECREF(args);
        if (!res)
            throw std::runtime_error("Call to ginput() failed.");

        const size_t len = PyList_Size(res);
        std::vector<std::array<double, 2>> out;
        out.reserve(len);
        for (size_t i = 0; i < len; i++) {
            PyObject* current = PyList_GetItem(res, i);
            std::array<double, 2> position;
            position[0] = PyFloat_AsDouble(PyTuple_GetItem(current, 0));
            position[1] = PyFloat_AsDouble(PyTuple_GetItem(current, 1));
            out.push_back(position);
        }
        Py_DECREF(res);

        return out;
    }

    inline void tight_layout() {
        PyObject* res = PyObject_CallObject(
            detail::_interpreter::get().s_python_function_tight_layout,
            detail::_interpreter::get().s_python_empty_tuple);

        if (!res)
            throw std::runtime_error("Call to tight_layout() failed.");

        Py_DECREF(res);
    }

    // recursion stop for the below
    template <typename... Args> bool plot() { return true; }

    // enable plotting of multiple triples (x, y, format)
    template <typename A, typename B, typename... Args>
    bool plot(const A& a, const B& b, const std::string& format, Args... args) {
        return plot(a, b, format) && plot(args...);
    }

    // This class allows dynamic plots, ie changing the plotted data without
    // clearing and re-plotting
    class Plot {
    public:
        // default initialization with plot label, some data and format
        template <typename Numeric>
        Plot(const std::string& name, const std::vector<Numeric>& x,
            const std::vector<Numeric>& y, const std::string& format = "") {

            assert(x.size() == y.size());

            PyObject* kwargs = PyDict_New();
            if (name != "")
                PyDict_SetItemString(kwargs, "label", PyString_FromString(name.c_str()));

            PyObject* xarray = get_array(x);
            PyObject* yarray = get_array(y);

            PyObject* pystring = PyString_FromString(format.c_str());

            PyObject* plot_args = PyTuple_New(3);
            PyTuple_SetItem(plot_args, 0, xarray);
            PyTuple_SetItem(plot_args, 1, yarray);
            PyTuple_SetItem(plot_args, 2, pystring);

            PyObject* res = PyObject_Call(
                detail::_interpreter::get().s_python_function_plot, plot_args, kwargs);

            Py_DECREF(kwargs);
            Py_DECREF(plot_args);

            if (res) {
                line = PyList_GetItem(res, 0);

                if (line)
                    set_data_fct = PyObject_GetAttrString(line, "set_data");
                else
                    Py_DECREF(line);
                Py_DECREF(res);
            }
        }

        // shorter initialization with name or format only
        // basically calls line, = plot([], [])
        Plot(const std::string& name = "", const std::string& format = "")
            : Plot(name, std::vector<double>(), std::vector<double>(), format) {}

        template <typename Numeric>
        bool update(const std::vector<Numeric>& x, const std::vector<Numeric>& y) {
            assert(x.size() == y.size());
            if (set_data_fct) {
                PyObject* xarray = get_array(x);
                PyObject* yarray = get_array(y);

                PyObject* plot_args = PyTuple_New(2);
                PyTuple_SetItem(plot_args, 0, xarray);
                PyTuple_SetItem(plot_args, 1, yarray);

                PyObject* res = PyObject_CallObject(set_data_fct, plot_args);
                if (res)
                    Py_DECREF(res);
                return res;
            }
            return false;
        }

        // clears the plot but keep it available
        bool clear() { return update(std::vector<double>(), std::vector<double>()); }

        // definitely remove this line
        void remove() {
            if (line) {
                auto remove_fct = PyObject_GetAttrString(line, "remove");
                PyObject* args = PyTuple_New(0);
                PyObject* res = PyObject_CallObject(remove_fct, args);
                if (res)
                    Py_DECREF(res);
            }
            decref();
        }

        ~Plot() { decref(); }

    private:
        void decref() {
            if (line)
                Py_DECREF(line);
            if (set_data_fct)
                Py_DECREF(set_data_fct);
        }

        PyObject* line = nullptr;
        PyObject* set_data_fct = nullptr;
    };

} // end namespace matplotlibcpp
