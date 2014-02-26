import sys
import matplotlib.pyplot as plt
import numpy as np

def plot_bounds(bounds, violated_constraints, output_name):
    # find when new lines start
    indices = list(zip(*bounds)[0])
    upper_bounds = list(zip(*bounds)[1])
    lower_bounds = list(zip(*bounds)[2])
    print("Found " + str(indices.count(0)) + " curves")
    curve_starts = list(np.where(np.array(indices) == 0)[0])
    curve_starts.append(len(indices) - 1)

    b = plt.figure()
    b.hold()

    for i in xrange(len(curve_starts) - 1):
    #for i in [0]:
        first_index = curve_starts[i]
        last_index = curve_starts[i + 1]

        # plot
        plt.plot(indices[first_index:last_index], upper_bounds[first_index:last_index], "*-")
        plt.plot(indices[first_index:last_index], lower_bounds[first_index:last_index], "^-")
        plt.legend(("Upper Bound", "Lower Bound"), loc='upper left', bbox_to_anchor=(1, 0.5))

        print("Plotting: ")
        print(indices[first_index:last_index])
        print(upper_bounds[first_index:last_index])
        print(lower_bounds[first_index:last_index])

    b.savefig("bounds_"+output_name+".png", dpi=800)

    c = plt.figure()
    for i in xrange(len(curve_starts) - 1):
    #for i in [0]:
        first_index = curve_starts[i]
        last_index = curve_starts[i + 1]

        # plot
        plt.plot(indices[first_index:last_index], violated_constraints[first_index:last_index], "o")
        plt.legend(("Violated Constraints"), loc='upper left', bbox_to_anchor=(1, 0.5))

    c.savefig("constraints_"+output_name+".png", dpi=800)


def parse_file(file, output_name):
    bounds = []
    violated_constraints = []

    for line in file:
        line = line.strip()

        if "Step" in line:
            # get the step number:
            step_number = line[line.rindex("Step  :") + 7:line.index("Value")].strip()

            # get the part between "Value :" and "Bound"
            value = line[line.rindex("Value :") + 7:line.index("Bound")].strip()

            # get the part after "Bound :"
            bound = line[line.rindex("Bound :") + 7:].strip()

            # store
            bounds.append([int(step_number), float(value), float(bound)])

        elif "Min number violated constraints: " in line:
            violated_constraints.append(int(line[line.rindex(" "):].strip()))

    plot_bounds(bounds, violated_constraints, output_name)




def main():
    print("Running tracking DD convergence plotter")
    if len(sys.argv) != 2:
        print("Wrong number or arguments given!\n")
        print("Usage: convergence_plot.py tracking_output.log")

        exit(0)

    filename = sys.argv[1]
    print("Reading file: " + filename)

    try:
        file = open(filename)
	output_name = filename.replace("output_","").replace(".log","")
        parse_file(file, output_name)

        # do something
    except IOError:
        print("Could not open file!")


if __name__ == "__main__":
    main()
