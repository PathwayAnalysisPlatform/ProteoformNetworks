import sys


def main():
    print(f"The number of arguments is: {len(sys.argv)}")
    print(f"Arguments: {sys.argv}")
    lines1 = set()
    lines2 = set()
    with open(sys.argv[1]) as f1:
        with open(sys.argv[2]) as f2:
            lines1 = {x.strip() for x in set(f1.readlines())}
            lines2 = {x.strip() for x in set(f2.readlines())}

    print(f"Len file set 1: {len(lines1)}")
    print(f"Len file set 2: {len(lines2)}")

    repeated = {x for x in lines1 if x in lines2}

    print(f"Num repeated: {len(repeated)}")
    print(f"Repeated: {repeated}")


if __name__ == "__main__":
    main()
