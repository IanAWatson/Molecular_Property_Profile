# Molecular Property Profile

A common need when examining collections of molecules is to be able
to compare the new molecules with one or more existing collections.
Answering questions like:

* How does the mean atom count differ between the collections
* Does this collection have more rotatable bonds
* How does `arbitrary molecular property` differ

The tool works in three phases.

1 Generate molecular features for all collections.
2 Use `generate_profile` to abstract the raw feature data into
 `Descriptor` proto form. This stores a summary of the data.
3. Use `plot_collections` to generate plots that superimpose
different distributions.

## Step 1
Any tool that can generate a tabular file of molecular properties
can be used. It must have a header. I have been using 
[iwdescr](https://github.com/IanAWatson/LillyMol-4.0-Bazel/blob/master/src/Molecule_Tools/iwdescr.cc)
which generates over 200 simple, and usually interpretable molecular
descriptors. But any such tool can be used. All that is needed is that
the result be tabular data, that can be read by Pandas.

In a large collection, computation may take a long time, the resulting
data file may be large, and these calculations may be awkward. Depending
on the need for completeness, large collections may be handled by
drawing random samples.

## Step 2
Convert these raw data files into protocol buffer form. The Descriptor
[proto](https...)
stores summary information about a feature, together with a description,
summary stats, and raw values.

If the number of distinct values is less than 100 (an arbitrarily chosen
number that can be changed), then all distinct values are stored. This
is the case even if the values are floats. If there are more than 100
distinct values, then an equally spaced sampling into 100 buckets is
done, using [numpy.histogram](https).

During generation of these protos, certain information must be specified.
The options currently supported are

* feature_name: the name(s) of the feature(s) being processed.
* feature_description: a meaningful description of what the feature is.
* collection: the source of these molecules (Pubchem, Chembl, you...)
* color: the color of any plots using this collection
* sep: the input file column separator.
* stem: a file name prefix that is used to generate proto files.

### Feature_name
This must be one of the column names in the input data file. The input is a
list, so enter multiple values in comma separated form: `natoms,nrings,aromdens`.

### feature_description
Optionally one can specify, for each feature, a description of that feature.
These are entered as separate options, so an invocation might look like
```
generate_profile --feature_name natoms,nrings,aromdens 
  --feature_description 'Number of heavy atoms' 
  --feature_description 'Number of rings' 
  --feature_description 'Aromatic density (fraction of atoms that are aromatic)'
```

If there are any feature_description options, there must be the same number as
the number of feature_name's specified. Or if this seems too tedious, the
resulting proto files can subsequently be edited.

### collection
A meaningful name for this souce of molecules. Commonly used values might include
Chembl, Pubchem, DrugBank, some molecules for sale, your virtual library,
your corporate collection...

### Color
Using a tool like this over many years showed the value of using consistent
colors for different collections - so that whenever people see a graph they know
that the green curve is Chembl. Select colors that seem meaningful and stick to
them.

### sep
This is passed to pd.read_csv. The default is space. Expect trouble if you
need to use things like tabs, in that case you might be better off changing
the default in the program.

### stem
Likely this will be the same as `collection`. The resulting proto files are
created with this as the leading part. You need to know which proto file(s)
came from Pubchem and which ones came from `your virtual library`.

For each feature_name specified, there will be a resulting `<stem>_<feature_name>.txt` file
created. For example, pubchem_natoms.txt. This contains a collection_pb2.Descrptor
proto describing the feature.

## Step 3
Once you have > 1 proto files generated, they can be graphically superimposed with
`plot_collections`. This is a command line tool that can either display
the plot on your terminal, or produces .png files if directed.

The individual collections previously profiled will all have different
X ranges. It is the job of plot_collections to line up those different
distributions so they can all be shown on the same scale.

If you have generated a bunch of distributions, you can automate
generating the .png files with something like

```
for feature in amw natoms htroaf ; do python plot_collections.py --stem STEM "CHEMBL_w_${feature}.dat" "RAND_w_${feature}.dat" "MUSH_w_${feature}.dat" ; done
```

