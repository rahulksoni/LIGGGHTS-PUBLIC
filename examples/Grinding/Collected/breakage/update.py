from paraview.simple import *
import glob

folder="post/dump*.liggghts"
stlfolder="post/dump*.stl"

def try_int(s):
        "Convert to integer if possible."
        try: return int(s)
        except: return s

def natsort_key(s):
        "Used internally to get a tuple by which s is sorted."
        import re
        return map(try_int, re.findall(r'(\d+|\D+)', s))

def natcmp(a, b):
        "Natural string comparison, case sensitive."
        return cmp(natsort_key(a), natsort_key(b))

def natcasecmp(a, b):
        "Natural string comparison, ignores case."
        return natcmp(a.lower(), b.lower())

# Find the source object
if (FindSource("liggghts_dump")==None): #create 1st time
	files = glob.glob(folder)
	files.sort(natcasecmp)
	xreader=liggghts_Reader(guiName="liggghts_dump", FileNames=files)
	#xreader.FileNames=files
	SetDisplayProperties(xreader, Representation="Point Sprite")
	DataRepresentation1 = Show()
	DataRepresentation1.PointSpriteDefaultsInitialized = 1
	DataRepresentation1.Texture = []
	DataRepresentation1.RadiusTransferFunctionEnabled = 1
	DataRepresentation1.RadiusMode = 'Scalar'
	DataRepresentation1.Representation = 'Point Sprite'
	DataRepresentation1.RadiusArray = [None, 'radius']
	DataRepresentation1.RadiusIsProportional = 1

else:
	files = glob.glob(folder) #update file list
	files.sort(natcasecmp)
	xreader.FileNames=files
	RenderView1 = GetRenderView()
	AnimationScene1 = GetAnimationScene()
	end=len(files)-1
	RenderView1.ViewTime = end
	AnimationScene1.AnimationTime = end
	Render()

# Find the source object
if (FindSource("stl_dump")==None): #create 1st time
	files = glob.glob(stlfolder)
	files.sort(natcasecmp)
	stlreader=STLReader(guiName="stl_dump", FileNames=files)
	
else:
	files = glob.glob(stlfolder) #update file list
	files.sort(natcasecmp)
	stlreader.FileNames=files
	RenderView1 = GetRenderView()
	AnimationScene1 = GetAnimationScene()
	end=len(files)-1
	RenderView1.ViewTime = end
	AnimationScene1.AnimationTime = end
	Render()
