import argparse
import matplotlib.pyplot as plt
import numpy as np
import os

directory = "curves"

def check(args):
	#
	if not os.path.exists(directory):
		os.makedirs(directory)
	#
	f_result = open(args.result, "r")
	f_answer = open(args.answer, "r")
	correct = True
	ans_arr = list()
	res_arr = list()
	#
	for idx, buff in enumerate(f_result):
		buff_ans = f_answer.readline()
		x, y = map(float, buff.split())
		res_arr.append(np.array([x, y]))
		a, b = map(float, buff_ans.split())
		ans_arr.append(np.array([a, b]))
		#
		chk_x = abs(x - a) <= 0.01
		chk_y = abs(y - b) <= 0.01
		#
		if not (chk_x and chk_y):
			print "Point #{0} is wrong!!!".format(idx)
			print " - Answer: x = {:.3f}, y = {:.3f}".format(a, b)
			print " - Result: x = {:.3f}, y = {:.3f}".format(x, y)
			correct = False
			#break
		count = idx
	#
	# Plot curve
	res_arr = np.vstack(res_arr)
	ans_arr = np.vstack(ans_arr)
	fig = plt.figure()
	res = plt.plot(res_arr[:, 0], res_arr[:, 1])
	fig.savefig("curves/{0}_res_curve.png".format(args.answer[:-4]))

	fig = plt.figure()
	ans = plt.plot(ans_arr[:, 0], ans_arr[:, 1])
	#plt.legend([ans, res], ("Answer", "Result"))
	fig.savefig("curves/{0}_ans_curve.png".format(args.answer[:-4]))
	#fig.savefig("curves/curve.png")
	#plt.show()
	#
	if correct:
		print "All case passed, {0} points checked".format(count + 1)
	f_result.close()
	f_answer.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("result", help="Result file.")
	parser.add_argument("answer", help="File to Compare.")
	
	args = parser.parse_args()
	
	check(args)


