

py:
	${MAKE} -C ./python

accel:
	${MAKE} -C ./accelpy

plot:
	python slabplot.py python/miniweather_mpi.slab

plotaccel:
	python slabplot.py accelpy/miniweather_accel.slab

tar:
	tar -cvf mini.tar *.png; rm -f *.png

profpy:
	${MAKE} -C ./python profile

profaccel:
	${MAKE} -C ./accelpy profile

clean:
	rm -f *.png *.tar
