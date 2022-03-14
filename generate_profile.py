"""Generates data for molecular property profiling.
"""

import sys

from absl import app
from absl import flags
from absl import logging

from google.protobuf import text_format

import numpy as np
import pandas as pd

from scipy import stats

import collection_pb2

# The max number of values that will be written
_MAX_FLOAT_POINTS = 100

FLAGS = flags.FLAGS

flags.DEFINE_list('feature_name', [], "Name of feature to process")
flags.DEFINE_multi_string('feature_description', [], 'Text description of each feature')
flags.DEFINE_string('collection', "", "name for the collection generating this data (Chembl...)")
flags.DEFINE_string('color', "black", "Line color when this collection is plotted")
flags.DEFINE_string('stem', "MPP", "name stem for files produced")
flags.DEFINE_boolean('verbose', False, "verbose output")

quantiles = [0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99]

def usage(ret):
  sys.exit(ret)

def determine_median(proto):
  """Set proto.median based on proto.quantile.

  Args:
  """
  proto.median = proto.quantile[3].value

def update_description(proto, collection: str,
                       feature_name:str,
                       feature_description: str,
                       line_color: str):
  """based on `name` update the description in `proto`.

  Args:
    proto:
    feature_name:
    feature_description:
  """
  proto.description.feature_name = feature_name

  if feature_description is not None:
    proto.description.plot_title = feature_description
  else:
    proto.description.plot_title = feature_name

  proto.description.line_color = line_color
  proto.description.source = collection
  if feature_description is not None:
    proto.description.description = feature_description

def add_quantiles(data: np.array,
                  proto: int):
  """
  """
  q = np.quantile(data, quantiles)
  for q, v in zip(quantiles, q):
    value = proto.quantile.add()
    value.quantile = q
    value.value = v

def set_numeric_values(data: np.array,
                       proto) -> None:
  """Update the numeric statistics in `proto` based on `data`.

  Args:
    data:
    proto:
  Returns:
  """
  proto.minval = float(np.min(data))
  proto.maxval = float(np.max(data))
  proto.mean = np.mean(data)
  add_quantiles(data, proto)
  determine_median(proto)

def profile_feature(data:np.array,
                    collection: str,
                    collection_color: str,
                    feature_name: str,
                    feature_description: str,
                    verbose: bool) -> int:
  """
  """
  if verbose:
    print(stats.describe(data))
  result = collection_pb2.Descriptor()
  update_description(result, collection, feature_name, feature_description, collection_color)

  set_numeric_values(data, result)

  unique, counts = np.unique(data, return_counts=True)
  if len(unique) < _MAX_FLOAT_POINTS:
    if data.dtype == np.int64:
      for v,c in zip(unique, counts):
        vc = result.int_values.add()
        vc.value = v
        vc.count = c
    else:
      for v,c in zip(unique, counts):
        vc = result.float_values.add()
        vc.value = v
        vc.count = c
    return result

  hist,bin_edges = np.histogram(data, bins=_MAX_FLOAT_POINTS)
  n = len(hist)
  for i in range(n):
#   if hist[i] == 0:
#     continue
    vc = result.float_values.add()
    vc.value = bin_edges[i]
    vc.count = hist[i]

  return result

def generate_feature_profile(data: pd.DataFrame,
                             collection: str,
                             feature_name: str,
                             feature_description: str,
                             collection_color: str,
                             name_stem: str,
                             verbose: bool) -> bool:
  """Generate a property profile for `feature_name`.


  Args:
    data:
    feature_name:
    feature_description:
    name_stem:
    verbose:
  Returns:
  """

  column_number = data.columns.get_loc(feature_name)
  if column_number < 0:
    logging.fatal("No %s in %r", feature_name, data.columns)

  feature_type = data.dtypes[column_number]

  if verbose:
    logging.info("Feature %s found in column %d type %r", feature_name, column_number, feature_type)

  proto = profile_feature(np.array(data[feature_name]), collection,
                         collection_color, feature_name, feature_description, verbose)

  output_fname = f"{name_stem}_{feature_name}.dat"
  with open(output_fname, "w") as writer:
    writer.write(text_format.MessageToString(proto))

def generate_profile(args):
  """Generates collection protos from molecular features.
  """
  verbose = FLAGS.verbose
  feature_name = FLAGS.feature_name
  feature_description = FLAGS.feature_description
  collection = FLAGS.collection
  collection_color = FLAGS.color
  name_stem = FLAGS.stem

  if len(args) == 1:
    logging.error("Must specify input file as argument")
    usage(1)

  if len(collection) == 0:
    logging.error("Must specify the collection name")
    usage(1)

  if len(feature_description) == 0:
    feature_description = [None] * len(feature_name)
  elif len(feature_name) != len(feature_description):
    logging.error("Must be a feature_description %d for each feature_name %d",
                   len(feature_description), len(feature_name))
    usage(1)

  if len(collection_color) == 0:
    collection_color = 'black'

  data = pd.read_csv(args[1], header=0, sep=' ')
  #data.set_index("Name", drop=True, inplace=True)
  if verbose:
    logging.info("Read dataframe with %d rows and %d columns", len(data), len(data.columns))

  if len(feature_name) == 0:
    feature_name = data.columns
  else:
    must_exit = False
    for name in feature_name:
      if not name in data.columns:
        logging.error("Cannot find %s in header", name)
        must_exit = True

    if must_exit:
      sys.exit(1)

  for i, name in enumerate(feature_name):
    generate_feature_profile(data, collection, name, feature_description[i],
                             collection_color, name_stem, verbose)



if __name__ == "__main__":
  app.run(generate_profile)
