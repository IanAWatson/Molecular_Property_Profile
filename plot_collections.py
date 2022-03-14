"""Plot one or more collections of molecule information
"""

from dataclasses import dataclass
import sys
from typing import List, Tuple

from absl import app
from absl import flags
from absl import logging

from google.protobuf import text_format

import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd

import collection_pb2

FLAGS = flags.FLAGS

flags.DEFINE_string("stem", "MPP", "name stem for .png files")
flags.DEFINE_float('qtruncate', 0.0, "Truncate plots at a quantile")
flags.DEFINE_boolean('verbose', False, "verbose output")

from dataclasses import dataclass

@dataclass
class Options:
    verbose: bool = False
    # The file name of plots generated
    stem: str = ""
    # We can truncate plots at a quantile
    truncate_quantile: float = 0.0

def usage(ret):
  sys.exit(ret)

def get_range(protos: List[collection_pb2.Descriptor]) -> Tuple:
  """Determine overall min and max of `protos`

  Args:
    protos: list of
  Returns:
  """

  # Get the range of the graph
  min_value = protos[0].minval
  max_value = protos[0].maxval
  for i in range(1, len(protos)):
    min_value = min(protos[i].minval, min_value)
    max_value = max(protos[i].maxval, max_value)
  logging.info("Range of %s btw %f and %f", protos[0].description.feature_name, min_value, max_value)
  return (min_value, max_value)

def int_plot(options: Options,
             protos: List,
             carray: List[np.array],
             varray: List[np.array]) -> None:
  """Generate a plot of integer data

  Args:
    options: options controlling behaviour
    protos: Descriptor protos
    carray: List of arrays of values extracted from `protos`
    varray: List of arrays of counts extracted from `protos`
  """
  (min_value, max_value) = get_range(protos)

  nvalues = round(max_value - min_value) + 1
  logging.info("Range %r to %r needs %d values", min_value, max_value, nvalues)

  # Each set will be aligned with this as the X axis
  x = np.arange(round(min_value), round(max_value) + 1)
  # A count array for each proto
  counts = []
  for proto in protos:
    c = np.zeros(nvalues)
    for vc in proto.int_values:
      ndx = round(vc.value - min_value)
      c[ndx] = vc.count
    counts.append(c)

  bars  = []
  width = 0.25
  for (i, c) in enumerate(counts):
    normed = c / np.linalg.norm(c, ord=1)
    b = plt.bar(x + width * i, normed, width,  color=protos[i].description.line_color)
    print(f"Setting label {protos[i].description.source}")
    b.set_label(protos[i].description.source)
    bars.append(b)

  plt.xlabel(protos[0].description.feature_name)
  plt.title(protos[0].description.description)
  plt.ylabel("Prevalence")
  plt.legend()
  plt.show()


def float_plot(options: Options,
               protos: List,
               carray: List[np.array],
               varray: List[np.array]) -> None:
  """Generate a plot of float data

  Args:
    options: options controlling behaviour
    protos: Descriptor protos
    carray: List of arrays of values extracted from `protos`
    varray: List of arrays of counts extracted from `protos`
  """
  (min_value, max_value) = get_range(protos)
  logging.info("Range %r to %r", min_value, max_value)

  # Each set will be aligned with this as the X axis
  dx = (max_value - min_value) / 50.0
  x = np.arange(min_value, max_value + dx, dx)
  nvalues = len(x)
  counts = []
  for proto in protos:
    c = np.zeros(nvalues)
    for vc in proto.float_values:
      ndx = round((vc.value - min_value) / dx)
      c[ndx] = vc.count
    counts.append(c)

  plts  = []
  width = 0.25
  for (i, c) in enumerate(counts):
    normed = c / np.linalg.norm(c, ord=1)
    p = plt.plot(x, normed, color=protos[i].description.line_color, label=protos[i].description.source)
    plts.append(p)

  plt.xlabel(protos[0].description.feature_name)
  plt.title(protos[0].description.description)
  plt.ylabel("Prevalence")
  plt.legend()
  plt.show()

def value_counts_to_arrays(from_proto) -> List[np.array]:
  """Convert the Int/Float ValueCount data `from_proto` to np arrays.

  Args:
    from_proto: source of data. Will be either IntValueCount or FloatValueCount
  Returns
    Two np arrays, one with the values, the other with counts.
  """
  n = len(from_proto)
  value = np.zeros(n)
  count = np.zeros(n)
  for i in range(n):
    vc = from_proto[i]
    value[i] = vc.value
    count[i] = vc.count

  return (value, count)

def do_plots(options: Options,
             protos:List) -> None:
  """Generate property profile plots from `protos`.

  Args:
    options:
    protos:
    name_stem:
  Returns:
  """
  varray:List[np.array] = []
  carray:List[np.array] = []
  is_int = 0
  for proto in protos:
    if len(proto.int_values) > 0:
      (value, count) = value_counts_to_arrays(proto.int_values)
      is_int += 1
    else:
      (value, count) = value_counts_to_arrays(proto.float_values)
    varray.append(value)
    carray.append(count)

  if is_int > 0:
    return int_plot(options, protos, varray, carray)
  else:
    return float_plot(options, protos, varray, carray)

def plot_profiles(args):
  """Generate plots for the protos in `args`.
  """
  options = Options()
  options.verbose = FLAGS.verbose
  options.stem = FLAGS.stem
  options.truncate_quantile = FLAGS.qtruncate

  if len(args) < 3:
    logging.error("Must specify at least two protos to plot")
    usage(1)

  protos = []
  feature_name = ""
  for i in range(1, len(args)):
    fname = args[i]
    with open(fname) as f:
      contents = f.readlines()

    if len(contents) == 0:
      logging.fatal("No data in %s", fname)

    proto = text_format.Parse('\n'.join(contents), collection_pb2.Descriptor())
    protos.append(proto)

    if feature_name == "":
      feature_name = proto.description.feature_name
    elif feature_name != proto.description.feature_name:
      logging.fatal("Name mismatch %s vs %s", feature_name, proto.feature_name)

  if options.verbose:
    logging.info("Read %d protos for property %s\n", len(protos), feature_name)

  do_plots(options, protos)


if __name__ == "__main__":
  app.run(plot_profiles)
