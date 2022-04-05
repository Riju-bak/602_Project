from batch import batchIVPSolver
from cont import continuosSolver


def displayBanner():
    print("--------------------------------------------------------")
    print("   WELCOME fellow Chemical Engineer :)   ")
    print(" Please choose a reactor type to simulate, or 3 to exit ")
    print("   1: Batch    2: Continuous    3: Exit   ")
    print("---------------------------------------------------------")


def main():
    while True:
        displayBanner()
        choice = input("Choice: ")
        if choice == '3':
            break
        elif choice == '1':
            batchIVPSolver()
        elif choice == '2':
            continuosSolver()
        else:
            print("Invalid Choice, input a number from 1,2, and 3.")


if __name__ == '__main__':
    main()
