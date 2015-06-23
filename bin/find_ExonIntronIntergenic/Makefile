cc=g++
exe=find_ExonIntronIntergenic
obj=main.o

$(exe):$(obj)
	$(cc) -o $(exe) $(obj)
main.o:main.cpp
	$(cc) -c main.cpp

clean:
	rm -rf *.o $(exe)