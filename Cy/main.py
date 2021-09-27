import pyximport; pyximport.install(reload_support=True, language_level=3)
import funcs
import time


t0 = time.time()
funcs.loop()
t1 = time.time()
print(t1 - t0)
