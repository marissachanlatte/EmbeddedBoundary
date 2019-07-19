from pathlib import Path
import numpy as np

def _split_and_remove_whitespace(line):
    return ' '.join(line.split()).split(' ')


def read_problem_specs(input_file):
    input_file = Path(input_file)
    assert input_file.exists(), "Input file: {} does not exist.".format(input_file)
    file = open(input_file, 'r')
    line = file.readline()
    size, start, stop, cells = _split_and_remove_whitespace(line)
    return float(size), float(start), float(stop), float(cells)


def build_mesh(size, start, stop, cells):
    x = np.linspace(start, stop, cells)
    y = np.linspace(start, stop, cells)
    X, Y = np.meshgrid(x, y)
    points = np.array(list(zip(X.flatten(), Y.flatten())))
    return points


def main():
    input_file = input("Enter path to input file: ")
    size, start, stop, cells = read_problem_specs(input_file)
    points = build_mesh(size, start, stop, cells)
    np.savetxt("mesh.txt", points)
    return 0

if __name__ == "__main__":
    main()
