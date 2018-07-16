with open("Makefile.am", "r") as fin:
    with open("makelist.txt", 'w') as fout:
        lines = fin.read().split('\n')
        for i, line in enumerate(lines):
            if line.find("_SOURCES")!=-1:
                print("Found")
                break
        for line in lines[i+1:]:
            if not len(line):
                break
            sl = line.split()[0:4]
            for l in sl:
                fout.write("<SourceFile>src/deps/clipper/{}</SourceFile>\n".format(l))
