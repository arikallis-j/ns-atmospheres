import sys, fire

from ns_atm.interface import *

def guide():
  print("This is guide.")
  
def main():
    if len(sys.argv) == 1:
        guide()
    else:
        calculator = Calculator()
        fire.Fire(calculator)

if __name__ == "__main__":
    main()
    