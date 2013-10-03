INCLUDE = -I/p/graphics/local/packages/fltk-1.0.11/include\
	 -I/p/graphics/local/packages/libtarga/include

LIB = -L/p/graphics/local/packages/fltk-1.0.11/lib\
	-L/p/graphics/local/packages/libtarga/lib\
	-L/usr/X11R6/lib

LINK = -lfltk -lX11 -lXext -ltarga

OBJ = ImageWidget.o ScriptHandler.o TargaImage.o 

Project1: $(OBJ)
	g++ -ggdb -Wall -o Project1 Main.cpp $(OBJ) $(INCLUDE) $(LIB) $(LINK) 

ImageWidget.o: ImageWidget.cpp ImageWidget.h
	g++ -ggdb -Wall -c -o ImageWidget.o ImageWidget.cpp $(INCLUDE)

ScriptHandler.o: ScriptHandler.cpp ScriptHandler.h
	g++ -ggdb -Wall -c -o ScriptHandler.o ScriptHandler.cpp $(INCLUDE)

TargaImage.o: TargaImage.cpp TargaImage.h
	g++ -ggdb -Wall -c -o TargaImage.o TargaImage.cpp $(INCLUDE)

clean:
	@for obj in $(OBJ); do\
		if test -f $$obj; then rm $$obj; fi; done
	@if (test -f Project1); then rm Project1; fi;

#libtarga.o: libtarga.c libtarga.h
#	gcc -c -o libtarga.o libtarga.c $(INCLUDE)
