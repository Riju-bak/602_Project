import sys
from batch import batchIVPSolver, generate_batch_yield_report, generate_batch_results_report
from cont import continuosSolver


def displayBanner():
    print("--------------------------------------------------------")
    print("   WELCOME fellow Chemical Engineer :)   ")
    print(" Please choose a reactor type to simulate, or 3 to exit ")
    print("   1: Batch    2: Continuous    3: Exit   ")
    print("---------------------------------------------------------")


def main(argv):
    if len(argv) == 2 and argv[1] == "gr":
        generate_batch_yield_report()
    elif len(argv) == 2 and argv[1] == "gbrr":
        generate_batch_results_report()
    elif len(argv) > 1:
        print("Usage: python main.py [gr]/[gbrr]")
    else:
        while True:
            displayBanner()
            choice = input("Choice: ")
            if choice == '3':
                break
            elif choice == '1':
                batchIVPSolver()
                generate_batch_results_report()
            elif choice == '2':
                continuosSolver()
            else:
                print("Invalid Choice, choose one out of 1,2, and 3.")


if __name__ == '__main__':
    main(sys.argv)
