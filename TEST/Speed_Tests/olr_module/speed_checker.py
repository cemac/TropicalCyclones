import time
import OLR

print('starting')
start_time = time.time()
c = OLR.OLR()
c.loop()
print("--- %s seconds ---" % (time.time() - start_time))
print('ending')
