default: all

all: release

debug:
	$(MAKE) -C debug

release:
	$(MAKE) -C release

profile:
	$(MAKE) -C profile

clean:
	$(MAKE) -C debug clean
	$(MAKE) -C release clean
	$(MAKE) -C profile clean

doc:
	doxygen Doxyfile
	
.PHONY: default all debug release profile clean doc