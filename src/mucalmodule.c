#include <Python.h>
#include <mucal.h>

static PyObject *MucalError;
static PyObject *mucal_mucal(PyObject *self, PyObject *args);

static PyMethodDef MucalMethods[] = {
  {"mucal", mucal_mucal, METH_VARARGS, "Calculated various x-ray parameters"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initmucal(void)
{
  PyObject *m;

  m = Py_InitModule("mucal", MucalMethods);
  if (m == NULL) return;

  MucalError = PyErr_NewException("mucal.error", NULL, NULL);
  Py_INCREF(MucalError);
  PyModule_AddObject(m, "error", MucalError);
}

static PyObject *
mucal_mucal(PyObject *self, PyObject *args)
{
  PyObject *ret;
  const char *element = "";
  int Z = 0;
  double energy = 0;

  double energies[9];
  double xsec[11];
  double fl_yield[4];

  char err_msg[100];
  int err;


  if (!PyArg_ParseTuple(args, "sd", &element, &energy) &&
      !PyArg_ParseTuple(args, "id", &Z, &energy))
  {
    return NULL;
  }

  err = mucal(
    element,
    Z,
    energy,
    'b',
    0,
    energies,
    xsec,
    fl_yield,
    err_msg
  );

  /*
  if (err)
  {
    PyErr_SetString(MucalError, err_msg);
    return NULL;
  }
  */

  ret = Py_BuildValue("{"
    "s:{s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d},"
    "s:{s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d},"
    "s:{s:d,s:d,s:d,s:d}"
    "}",
    "energies",
    "K",
    energies[0],
    "L1",
    energies[1],
    "L2",
    energies[2],
    "L3",
    energies[3],
    "M",
    energies[4],
    "Ka",
    energies[5],
    "Kb",
    energies[6],
    "La",
    energies[7],
    "Lb",
    energies[8],
    "xsec",
    "photoelectric",
    xsec[0],
    "coherent",
    xsec[1],
    "incoherent",
    xsec[2],
    "total",
    xsec[3],
    "conversion",
    xsec[4],
    "abscoeff",
    xsec[5],
    "atwt",
    xsec[6],
    "density",
    xsec[7],
    "l1_jump",
    xsec[8],
    "l2_jump",
    xsec[9],
    "l3_jump",
    xsec[10],
    "fl_yield",
    "K",
    fl_yield[0],
    "L1",
    fl_yield[1],
    "L2",
    fl_yield[2],
    "L3",
    fl_yield[3]
    );

  return ret;
}
