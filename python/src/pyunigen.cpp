// Python bindings for Unigen, heavily based on the Python bindings written for CryptoMiniSat.

#include <Python.h>
#include <unigen/unigen.h>
#include <approxmc/approxmc.h>
#include <cryptominisat5/cryptominisat.h>

#include <limits>

#define MODULE_NAME "pyunigen"
#define MODULE_DOC "Unigen almost uniform sampler."

typedef struct {
    PyObject_HEAD
    UniGen::UniG* unig;
    ApproxMC::AppMC* appmc;
    std::vector<CMSat::Lit> tmp_cl_lits;
    PyObject* sampled_val;

    int verbosity;
    uint32_t seed;
    double kappa;
    std::vector<uint32_t> sampling_set;
} Sampler;

static const char sampler_create_docstring[] = \
"Sampler(verbosity=0, seed=1, kappa=0.638, sampling_set=None)\n\
Create Sampler object.\n\
\n\
:param verbosity: Verbosity level: 0: nothing printed; 15: very verbose.\n\
:param seed: Random seed\n\
:param kappa: Uniformity parameter (see TACAS-15 paper)\n\
:param sampling_set: (Optional) If provided, the solutions are sampled almost uniformly\n\
    over the variables in sampling_set.";

/********** Internal Functions **********/

/* Helper functions */

void pybinding_callback(const std::vector<int>& solution, void *source)
{
    PyObject** py_obj_cont = (PyObject**) source;

    PyObject* py_list = PyList_New(solution.size());
    if (py_list == NULL) {
        PyErr_SetString(PyExc_SystemError, "failed to create a list");
        return;
    }
    for (unsigned int i = 0; i < solution.size(); i++) {
        PyObject *lit = PyLong_FromLong((long)solution[i]);
        if (lit == NULL) {
            PyErr_SetString(PyExc_SystemError, "failed to create a list");
            return;
        }
        PyList_SET_ITEM(py_list, i, lit);
    }

    (*py_obj_cont) = py_list;
}

static int parse_sampling_set(Sampler *self, PyObject *sample_set_obj)
{
    PyObject *iterator = PyObject_GetIter(sample_set_obj);
    if (iterator == NULL) {
        PyErr_SetString(PyExc_TypeError, "iterable object expected");
        return 1;
    }

    PyObject *lit;
    while ((lit = PyIter_Next(iterator)) != NULL) {
        long val = PyLong_AsLong(lit);
        if (val == 0) {
            PyErr_SetString(PyExc_ValueError, "non-zero integer expected");
            return 1;
        }
        if (val > std::numeric_limits<int>::max()/2
            || val < std::numeric_limits<int>::min()/2
        ) {
            PyErr_Format(PyExc_ValueError, "integer %ld is too small or too large", val);
            return 1;
        }

        long var = val - 1;
        self->sampling_set.push_back(var);
        Py_DECREF(lit);
    }
    Py_DECREF(iterator);

    return 0;
}

static void setup_sampler(Sampler *self, PyObject *args, PyObject *kwds)
{
    static char* kwlist[] = {"verbosity", "seed", "kappa", "sampling_set", NULL};

    // All parameters have the same default as the command line defaults
    // except for verbosity which is 0 by default.
    self->verbosity = 0;
    self->seed = 1;
    self->kappa = 0.638;

    PyObject* sample_set_obj = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iIdO", kwlist,
        &self->verbosity, &self->seed, &self->kappa, &sample_set_obj))
    {
        return;
    }

    if (sample_set_obj != NULL && parse_sampling_set(self, sample_set_obj)) {
        return;
    }

    if (self->verbosity < 0) {
        PyErr_SetString(PyExc_ValueError, "verbosity must be at least 0");
        return;
    }
    if (self->kappa <= 0 || self->kappa >= 1) {
        PyErr_SetString(PyExc_ValueError, "epsilon must be greater than 0");
        return;
    }

    self->appmc = new ApproxMC::AppMC;
    self->appmc->set_verbosity(self->verbosity);
    self->appmc->set_seed(self->seed);
    self->appmc->set_projection_set(self->sampling_set);

    self->unig = new UniGen::UniG(self->appmc);
    self->unig->set_verbosity(self->verbosity);
    self->unig->set_kappa(self->kappa);
    // Only one sample generated and
    // multisample disabled for simplicity
    self->unig->set_multisample(false);

    self->unig->set_callback(pybinding_callback, &self->sampled_val);

    return;
}

static int convert_lit_to_sign_and_var(PyObject* lit, long& var, bool& sign)
{
    if (!PyLong_Check(lit))  {
        PyErr_SetString(PyExc_TypeError, "integer expected !");
        return 0;
    }

    long val = PyLong_AsLong(lit);
    if (val == 0) {
        PyErr_SetString(PyExc_ValueError, "non-zero integer expected");
        return 0;
    }
    if (val > std::numeric_limits<int>::max()/2
        || val < std::numeric_limits<int>::min()/2
    ) {
        PyErr_Format(PyExc_ValueError, "integer %ld is too small or too large", val);
        return 0;
    }

    sign = (val < 0);
    var = std::abs(val) - 1;

    return 1;
}

static int parse_clause(Sampler *self, PyObject *clause, std::vector<CMSat::Lit>& lits)
{
    PyObject *iterator = PyObject_GetIter(clause);
    if (iterator == NULL) {
        PyErr_SetString(PyExc_TypeError, "iterable object expected");
        return 0;
    }

    PyObject *lit;
    long int max_var = 0;
    while ((lit = PyIter_Next(iterator)) != NULL) {
        long var;
        bool sign;
        int ret = convert_lit_to_sign_and_var(lit, var, sign);
        Py_DECREF(lit);
        if (!ret) {
            Py_DECREF(iterator);
            return 0;
        }
        max_var = std::max(var, max_var);

        lits.push_back(CMSat::Lit(var, sign));
    }

    if (!lits.empty() && max_var >= (long int)self->appmc->nVars()) {
        self->appmc->new_vars(max_var-(long int)self->appmc->nVars()+1);
    }

    Py_DECREF(iterator);
    if (PyErr_Occurred()) {
        return 0;
    }

    return 1;
}

static int _add_clause(Sampler *self, PyObject *clause)
{
    self->tmp_cl_lits.clear();
    if (!parse_clause(self, clause, self->tmp_cl_lits)) {
        return 0;
    }
    self->appmc->add_clause(self->tmp_cl_lits);

    return 1;
}

/* add_clause function */

PyDoc_STRVAR(add_clause_doc,
"add_clause(clause)\n\
Add a clause to the solver.\n\
\n\
:param clause: An iterator contains literals (ints)"
);

static PyObject* add_clause(Sampler *self, PyObject *args, PyObject *kwds)
{
    static char* kwlist[] = {"clause", NULL};
    PyObject *clause;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &clause)) {
        return NULL;
    }

    if (_add_clause(self, clause) == 0 ) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;

}

/* sample function */

PyDoc_STRVAR(sample_doc,
"sample()\n\
Sample almost uniformly from the solutions for the clauses that have been \n\
added with add_clause().\n\
\n\
:param cell_count: The number of solutions in a cell.\n\
:param hash_count: The number of hashes applied to the formula.\n\
:return: A Python list containing a sampled solution."
);

static PyObject* sample(Sampler *self, PyObject *args, PyObject *kwds)
{
    self->appmc->setup_vars();

    self->sampled_val = NULL;

    static char* kwlist[] = {"cell_count", "hash_count", NULL};

    auto sol_count = * new ApproxMC::SolCount;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "II", kwlist,
        &sol_count.cellSolCount, &sol_count.hashCount))
    {
        return NULL;
    }

    self->unig->sample(&sol_count, 1);

    if (self->sampled_val == NULL) {
        return NULL;
    }

    return self->sampled_val;
}

/********** Python Bindings **********/
static PyMethodDef Sampler_methods[] = {
    {"sample",    (PyCFunction) sample,      METH_VARARGS | METH_KEYWORDS, sample_doc},
    {"add_clause",(PyCFunction) add_clause,  METH_VARARGS | METH_KEYWORDS, add_clause_doc},
    {NULL, NULL}  // Sentinel
};

static void Sampler_dealloc(Sampler* self)
{
    delete self->unig;
    delete self->appmc;
    Py_TYPE(self)->tp_free ((PyObject*) self);
}

static int Sampler_init(Sampler *self, PyObject *args, PyObject *kwds)
{
    if (self->unig != NULL) {
        delete self->unig;
    }

    if (self->appmc != NULL) {
        delete self->appmc;
    }

    setup_sampler(self, args, kwds);

    if (!self->unig) {
        return -1;
    }

    if (!self->appmc) {
        return -1;
    }

    return 0;
}

static PyTypeObject pyunigen_SamplerType = 
{
    PyVarObject_HEAD_INIT(NULL, 0)  /*ob_size*/
    "pyapproxmc.Sampler",           /*tp_name*/
    sizeof(Sampler),                /*tp_basicsize*/
    0,                              /*tp_itemsize*/
    (destructor)Sampler_dealloc,    /*tp_dealloc*/
    0,                              /*tp_print*/
    0,                              /*tp_getattr*/
    0,                              /*tp_setattr*/
    0,                              /*tp_compare*/
    0,                              /*tp_repr*/
    0,                              /*tp_as_number*/
    0,                              /*tp_as_sequence*/
    0,                              /*tp_as_mapping*/
    0,                              /*tp_hash */
    0,                              /*tp_call*/
    0,                              /*tp_str*/
    0,                              /*tp_getattro*/
    0,                              /*tp_setattro*/
    0,                              /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    sampler_create_docstring,       /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    Sampler_methods,                /* tp_methods */
    0,                              /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    (initproc)Sampler_init,         /* tp_init */
};

PyMODINIT_FUNC PyInit_pyunigen(void) 
{
    PyObject* m;

    pyunigen_SamplerType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&pyunigen_SamplerType) < 0) {
        // Return NULL on Python3 and on Python2 with MODULE_INIT_FUNC macro
        // In pure Python2: return nothing.
        return NULL;
    }

    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,  /* m_base */
        MODULE_NAME,            /* m_name */
        MODULE_DOC,             /* m_doc */
        -1,                     /* m_size */
        NULL,                   /* m_methods */
        NULL,                   /* m_reload */
        NULL,                   /* m_traverse */
        NULL,                   /* m_clear */
        NULL,                   /* m_free */
    };

    m = PyModule_Create(&moduledef);

    if (!m) {
        return NULL;
    }

    // Add the Sampler type
    Py_INCREF(&pyunigen_SamplerType);
    if (PyModule_AddObject(m, "Sampler", (PyObject *)&pyunigen_SamplerType)) {
        Py_DECREF(m);
        return NULL;
    }

    return m;
}

