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

If you decide to use iwdescr, the correct invocation would be
```
iwdescr -l -g all -O all -i ICTE -i smi file.smi > file.dat
```

The `-i ICTE` tells it to ignore connection table errors. LillyMol is
not very good at finding Kekule forms for aromatic smiles, and so may
fail, which would otherwise terminate the computation.

## Step 2
Convert these raw data files into protocol buffer form. The Descriptor
[proto](https...)
stores summary information about a feature, together with a description,
summary stats, and sampled values.

If the number of distinct values is less than 100 (an arbitrarily chosen
number that can be changed), then all distinct values are stored. This
is the case even if the values are floats. If there are more than 100
distinct values, then an equally spaced sampling into 100 buckets is
done, using [numpy.histogram](https).

For each feature to be plotted, it is useful to have a description
of what that feature is. So while a column name might be `natoms` a
more meaningful name for a plot might be `Number of Heavy Atoms`. I 
have included a file, `column_descriptions.txt`, which contains
column name descriptions for some of the columns generated by `iwdescr`.

Given that file, the easy way to then generate a summary proto
for calculated features in a file is via something like
(for a collection called `rand`):

```
#python generate_profile.py --read_descriptions=column_descriptions.txt
                            --collection rand
                            --stem RAND
                            --color=red
                            rand.dat
```

This processes all columns in `column_descriptions.txt` and will generate
summary data for each column. Each feature is written to a separate file
of the form `<stem>_<feature_name>.txt'

Precomputed files for Chembl, Pubchem and Enamine HTS are included.

Options:

### read_descriptions
The name of an existing file that contains feature descriptions.
During generation of these protos, certain information must be specified.
The `column_descriptions.txt` can be used.

### collection
The name of the collection, likely something like Pubchem, Chembl,
molecules from vendor, my virtual library, company collection ....

### stem
A file name stem for the files produced. This should likely be the
same as the collection name.

### color
Using a tool like this over many years showed the value of using consistent
colors for different collections - so that whenever people see a graph they know
that the green curve is Chembl. Select colors that seem meaningful and stick to
them.

### sep
This is passed to pd.read_csv. The default is space. Expect trouble if you
need to use things like tabs, in that case you might be better off changing
the default in the program.

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
for feature in amw natoms htroaf ; do python plot_collections.py --stem STEM "CHEMBL_w_${feature}.dat" "RAND_w_${feature}.dat" "PUBCHEM_w_${feature}.dat" ; done
```

Or some clever wildcarding can be used to process all features.

The tool supports some options

### stem
Generate a .png file of the form `<stem>_<feature>.png` rather than displaying
on the screen.

### X and Y
Specify a plot size in inches, `--X 6.0 --Y 8.0`. Either omit both, or they
must both be specified.

### xmin and xmax
Specify the x axis range. By default, it will be such that all values from
all collections are displayed. So, if there were a collection that had some
very large molecules, it might be clearer to see the relevant molecules
via something like `-xmin 0 -xmax 60` when viewing heavy atom distributions.

# Summary
A complete workflow, for the collection `foo` might look like
```

iwdescr -i ICTE -i smi -g all -l -O all foo.smi > foo.dat

# Check for errors

python generate_profile.py --read_descriptions=column_descriptions.txt
                           --collection foo
                           --stem FOO
                           --color=red
                            foo.dat

# Check to make sure all FOO*.dat files are produced

# Generate comparison plots for rotatable bonds.

python plot_collections.py FOO_w_rotbond.dat CHEMBL_w_rotbond.dat ENAMINE_w_rotbond.dat

# Or generate plots for all features, based on what is in `column_descriptions.txt`

plot_collections.py --quantile 0.01 --collection FOO,CHEMBL,ENAMINE --feature_description column_descriptions.txt

# Note the quantile argument to remove extreme outliers, which can
# make interpretation more difficult.
```

This will display the plots on the terminal. Add a `--stem` argument and instead
plots will be written to .png files.
