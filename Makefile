default: all

all: debug release profile

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
	
.PHONY: default all debug release profile clean