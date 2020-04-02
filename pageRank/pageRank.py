import os
import sys

def ceil(a):
	if int(a) == a:
		return int(a)
	else:
		return int(a) + 1


def deepCopy(dic_src, dic_dst):
	dic_dst.clear()
	for key in dic_src:
		dic_dst[key] = dic_src[key]


def makeBlock(path, indexList):
	sortedIndexList = sorted(indexList)
	blockSize = min(1000, len(indexList))
	
	blockNum = ceil(len(sortedIndexList) / blockSize)
	
	upperBound = sortedIndexList[blockSize - 1]
	lowerBound = sortedIndexList[0]
	blockCount = 0
	while True:
		with open(path, 'r') as f:
			inversedDic = dict()
			base = blockCount * blockSize
			for num in sortedIndexList[base: base + blockSize]:
				inversedDic[num] = []
			string = f.readline()
			while string != '':
				data = [int(x) for x in string.split('\t')]
				if data[0] >= lowerBound and data[0] <= upperBound:
					inversedDic[data[0]].append(data[1])
				string = f.readline()
			storeBlock("_record/block_" + str(blockCount), inversedDic)
		blockCount += 1
		base = blockSize * blockCount
		if blockCount == blockNum:
			break
		elif blockCount == blockNum - 1:
			upperBound = sortedIndexList[len(sortedIndexList) - 1]
		else:
			upperBound = sortedIndexList[base + blockSize - 1]
		lowerBound = sortedIndexList[base]


def storeBlock(name, dic):
	with open(name + ".txt", 'w') as f:
		for key in dic:
			val = dic[key]
			val.insert(0, key)
			line = ' '.join([str(x) for x in val])
			f.write(line + '\n')


def loadBlock(name):
	dic = dict()
	with open(name + '.txt', 'r') as f:
		lines = f.readlines()
		for line in lines:
			data = line.strip('\n').split(' ')
			data = [int(x) for x in data]
			if len(data) > 1:
				dic[data[0]] = data[1:]
			else:
				dic[data[0]] = []
	return dic


def makeIndex(path):
	l = []
	with open(path, 'r') as f:
		str = f.readline()
		while str != '':
			data = [int(x) for x in str.strip('\n').split('\t')]
			if data[0] not in l:
				l.append(data[0])
			if data[1] not in l:
				l.append(data[1])
			str = f.readline()
	return l


def storeIndex(name, l):
	with open(name + ".txt", 'w') as f:
		f.write(' '.join([str(x) for x in sorted(l)]))


def loadIndex(name):
	with open(name + '.txt', 'r') as f:
		line = f.readline()
		return [int(x) for x in line.split(' ')]


def preprocess(path):
	if not os.path.exists('_record'):
		os.mkdir("_record")
		print("Initial record directory \"_record\"\n")
	else:
		print("Record directory \"_record\" already exists\n")
	print("Building index...\n")
	index = makeIndex(path)
	storeIndex('_record/index', index)
	print("Index building done\n")
	print("Dividing data into block...\n")
	makeBlock(path, index)
	print("Data reconstruction done\n")


def deleteMediumRecord():
	if os.path.exists("_record"):
		for file in os.listdir("_record"):
			os.remove("_record/" + file)
		os.removedirs("_record")


def parseArgs():
	args = len(sys.argv)
	argv = sys.argv
	if args == 1:
		print("you need to assign source path")
		print("Run with parameter \"help\" to get more info")
		sys.exit(-1)
	elif args == 2:
		if argv[1] == "help":
			print("one mandatory parameter and three optional parameters are needed")
			print("source path(madatory)")
			print("beta(optional): random teleport possibility")
			print("epsilon(optional): exponent of convergence condition")
			print("delete medium record(optional): T/F")
			sys.exit(0)
		else:
			if os.path.exists(argv[1]):
				print("Run with default parameter")
				return [argv[1], 0.8, 10, True]
			else:
				print("source path \"" + argv[1] + "\" doesn't exist")
				sys.exit(-1)
	elif args == 5:
		if os.path.exists(argv[1]):
			b = 0.8
			e = 10
			isDelete = True
			try:
				b = float(argv[2])
			except:
				print("Invaild parameter \"" + argv[2] + "\"")
				sys.exit(-1)
			try:
				e = float(argv[3])
			except:
				print("Invaild parameter \"" + argv[3] + "\"")
				sys.exit(-1)
			if argv[4] == 'T':
				isDelete = True
			elif argv[4] == 'F':
				isDelete = False
			else:
				print("Invaild parameter \"" + argv[3] + "\"")
				sys.exit(-1)
			return [argv[1], b, e, isDelete]
		else:
			print("source path \"" + argv[0] + "\" doesn't exist")
			sys.exit(-1)
	else:
		print("Incorrect parameter number")
		sys.exit(-1)


def norm2dis(dic1, dic2):
	distance = 0
	for key in dic1.keys():
		distance += (dic1[key] - dic2[key])**2
	return distance


def normalization(dic):
	sum = 0.0
	for key in dic:
		sum += dic[key]
	for key in dic:
		dic[key] /= sum


def pageRank(b, e, recordRemove):
	print("Start calculating page rank score...\n")
	print("parameters: \nbeta\t" +str(b)+ "\nepsilon\t10^(-" +str(e)+ ")")
	print("Delete medium Record\t" +str(recordRemove)+ "\n")
	beta = b
	epsilon = 10 ** (-e)
	index = loadIndex('_record/index')
	size = len(index)
	modifyFactor = 1 / size
	blockSize = min(1000, size)
	blockNum = ceil(size / blockSize)
	iter = 0
	
	currDic = dict()
	oldDic = dict()
	

	for i in index:
		currDic[i] = 1 / size
		oldDic[i] = 1
	
	while norm2dis(currDic, oldDic) > epsilon:
		iter += 1
		deepCopy(currDic, oldDic)
		blockCount = 0
		while blockCount < blockNum:
			block = loadBlock('_record/block_'+str(blockCount))
			for key in block:
				dst = block[key]
				if dst != []:
					updateVal = beta * (oldDic[key] / len(dst)) \
					            + (1 - beta) * modifyFactor
					for d in dst:
						currDic[d] += updateVal
			blockCount += 1
		normalization(currDic)
	with open('pageRank.txt', 'w') as f:
		lines = []
		for key in currDic:
			line = str(key) + ': ' + str(currDic[key]) + '\n'
			lines.append(line)
		f.writelines(lines)
	if recordRemove:
		deleteMediumRecord()
	print("Page Rank done\nTotal iteration times: " + str(iter) + "\n")

if __name__ == "__main__":
	argv = parseArgs()
	preprocess(argv[0])
	pageRank(argv[1], argv[2], argv[3])