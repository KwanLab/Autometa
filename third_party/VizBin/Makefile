
JAR = dist/VizBin-dist.jar
.PHONY: VizBin clean

VizBin:
	ant

install: $(JAR)
install:
	@echo "COPYING: jar into root folder ..."
	install -m 0755 $(JAR) ../../../
	@echo "done."

clean:
	@echo "Nothing to be cleaned in VizBin."
