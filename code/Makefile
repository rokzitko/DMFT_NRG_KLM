.PHONY: run clean

run: DONE

clean:
	scripts/cleanup
	rm -rf __pycache__ scripts/__pycache__

DONE: template/data.in DOS.dat PHI.dat
	./START
