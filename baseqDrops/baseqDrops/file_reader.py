import subprocess

def read_file_by_lines(filepath, maxLines, linecount):
    counter = 0
    if filepath.endswith("gz") or filepath.endswith("gzip") or filepath.endswith("gz2"):
        reading = subprocess.Popen(["gunzip", "-c", filepath], stdout=subprocess.PIPE, bufsize=1000000)
        infile = reading.stdout
        while True:
            counter = counter + 1
            if counter > maxLines:
                return
            data = [infile.readline().decode('utf8') for i in range(linecount)]
            if data[0] == "":
                return
            yield data
    else:
        infile = open(filepath, 'r')
        while True:
            counter = counter + 1
            if counter > maxLines:
                return
            data = [infile.readline() for i in range(linecount)]
            if data[0] == "":
                return
            yield data

def read_filelines(filepath, maxLines, linecount, skip=0):
    counter = 0
    maxLines = maxLines+skip
    if filepath.endswith("gz") or filepath.endswith("gzip") or filepath.endswith("gz2"):
        reading = subprocess.Popen(["gunzip", "-c", filepath], stdout=subprocess.PIPE, bufsize=1000000)
        infile = reading.stdout
        while True:
            counter = counter + 1
            if counter <= skip:
                continue
            if counter > maxLines:
                return
            data = [infile.readline().decode('utf8') for i in range(linecount)]
            if data[0] == "":
                return
            yield data
    else:
        infile = open(filepath, 'r')
        while True:
            counter = counter + 1
            if counter <= skip:
                continue
            if counter > maxLines:
                return
            data = [infile.readline() for i in range(linecount)]
            if data[0] == "":
                return
            yield data