"""Executes the main steps of conversion"""
import os
import sys
import warnings
import argparse
import read_njoy_outputs
import combiner
import xs_writer


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.MetavarTypeHelpFormatter,
                      argparse.RawTextHelpFormatter):
    pass


# ------------------------------------------------------------
# Check Python version
# ------------------------------------------------------------

if sys.version_info[0] < 3:
    print("\n Error: This script requires python3 but was executed "
          "with version:\n\n"+sys.version+"\n")
    sys.exit(1)

# ------------------------------------------------------------
# Create parser
# ------------------------------------------------------------

argparser = argparse.ArgumentParser(
    description="Executes the main steps of conversion",
    formatter_class=CustomFormatter
)

argparser.add_argument(
    '--output-directory',
    type=str,
    required=True,
    help="Complete path where to store the output.",
)

argparser.add_argument(
    '--njoy-output-filename',
    type=str,
    required=True,
    help="Name of output file produced by NJOY."
)

argparser.add_argument(
    '--xs-filename',
    type=str,
    help="Name of XS file.",
)

argparser.add_argument(
    '--plot',
    action='store_true',
    default=False,
    help="A flag for plotting the XS data."
)

argv = argparser.parse_args()

output_directory = os.path.abspath(argv.output_directory)
if not os.path.isdir(output_directory):
    warnings.warn(
        "Value supplied to --output-directory does not point to an "
        f"existing directory. Creating one at {output_directory}."
    )
    os.makedirs(output_directory)

# ------------------------------------------------------------
# Read NJOY output
# ------------------------------------------------------------

njoy_output_path = os.path.join(
    output_directory, "njoy", argv.njoy_output_filename
)
print(f"\nReading NJOY output located at {njoy_output_path}...")

raw_njoy_data = read_njoy_outputs.read_njoy_file(
    njoy_output_path, verbose=True
)

# ------------------------------------------------------------
# Combine raw data
# ------------------------------------------------------------

print("\nCombining raw data...")

data, problem_description = combiner.build_combined_data(
    raw_njoy_data, plot=argv.plot
)

# TODO: this is not a good way of recovering the isotope if the user
#       supplies its own output filename
filename = argv.njoy_output_filename.split('_')
problem_description['isotope'] = filename[0].split('.')[0]
print(f"\nProblem Description:")
for key, val in problem_description.items():
    print(f"{key:<15}: {val}")

# ------------------------------------------------------------
# Write cross-section file
# ------------------------------------------------------------

if len(argv.xs_filename) == 0:
    # generate name automatically
    xs_filename = f"{argv.njoy_output_filename}.xs"
else:
    xs_filename = argv.xs_filename

xs_output_directory = os.path.join(output_directory, "xs")
if not os.path.isdir(xs_output_directory):
    warnings.warn(f"Output directory for `.xs` files does not exist. "
                  f"Creating it at {xs_output_directory}.")
    os.makedirs(xs_output_directory)

xs_output_path = os.path.join(xs_output_directory, xs_filename)
print(f"\nCreating cross-section in file {xs_filename}")

xs_writer.write_xs_file(
    data, xs_output_path, problem_description
)

# ------------------------------------------------------------
# Move plots
# ------------------------------------------------------------

if argv.plot:
    import os
    import glob

    filename = xs_filename.split('.')
    for file in glob.glob("*.png"):
        plot_dir = os.path.join(output_directory, "figures")
        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir)
        figure_path = os.path.join(plot_dir, f"{filename[0]}_{file}")
        os.system(f"mv {file} {figure_path}")
