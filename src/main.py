from batch import batchSolver, batchIVPSolver, generate_yield_report, generate_batch_results_report
import sys, getopt


def main(argv):
    if len(argv) == 2 and argv[1] == "gr":
        generate_yield_report()
    elif len(argv)==2 and argv[1] == "gbrr":
        generate_batch_results_report()
    elif len(argv) > 1:
        print("Usage: python main.py [gr]/[gbrr]")
    else:
        batchIVPSolver()


if __name__ == '__main__':
    main(sys.argv)