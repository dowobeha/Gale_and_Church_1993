
all: file1.al file2.al
	cat file1.al file2.al

file1.al file2.al: align_regions file1 file2
	./align_regions -D '.PARA' -d '.End of Sentence' -v file1 file2

align_regions: align_regions.c
	cc -std=c89 align_regions.c -lm -o align_regions 

clean:
	rm -f file1.al file1.al align_regions

