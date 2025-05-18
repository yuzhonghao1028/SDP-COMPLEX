import os  
import sys  
import json  
import time  
import numpy as np  
import vtk  
import tkinter as tk  
from tkinter import filedialog  
from tqdm import tqdm  
from shapely.geometry import Point  # 引入Shapely库  

# ───────────────────── ① 输入用户文件夹 ─────────────────────  
def get_user_input():  
    root = tk.Tk()  
    root.withdraw()  # 隐藏主窗口  

    base_dir = filedialog.askdirectory(title="请选择文件夹")  
    if not base_dir:  
        sys.exit("[Err] 文件夹未选择！")  

    return base_dir  

# ───────────────────── ② 固定参数 ─────────────────────  
EXTEND_LEN = 20.0  # 起点反向外延设置为20mm  
SAMPLE_STEP = 0.01  # 固定步长为0.01 mm  
PLANE_SKIP = 1  # 设置取点的跳过数为1（即不跳过任何点）  

IGNORE_R = 3.0  
TUBE_DIAM = 3.0  

# ───────────────────── ③ 固定目录 ─────────────────────  
BASE_DIR = get_user_input()  # 获取用户选择的目录  

entry_path  = os.path.join(BASE_DIR, "entry.vtk")  
hard_path   = os.path.join(BASE_DIR, "nonSurgicalSpace.vtk")  

ENDPOINTS = [  
    ("anterior.mrk.json",  "anterior",  (1.0, 1.0, 0.0)),   
    ("middle.mrk.json",    "middle",    (13/255, 1.0, 0.0)),  
    ("posterior.mrk.json", "superior",  (1.0, 0.0, 0.0)),   
]  

# ───────────────────── ④ 工具函数 ─────────────────────  
def load_poly(path):  
    if not os.path.exists(path):  
        sys.exit(f"[Err] 文件不存在: {path}")  
    rd = vtk.vtkPolyDataReader()  
    rd.SetFileName(path)  
    rd.Update()  
    pd = rd.GetOutput()  
    if pd.GetNumberOfPoints() == 0:  
        sys.exit(f"[Err] 模型为空: {path}")  
    return pd  

def poly_points(pd):  
    g = pd.GetPoints()  
    return [g.GetPoint(i) for i in range(g.GetNumberOfPoints())]  

def read_endpoints(js_path):  
    cps = json.load(open(js_path, encoding="utf-8"))['markups'][0]['controlPoints']  
    return [cp['position'] for cp in cps]  

def hits_cylinder(p0, p1, tree, radius):  
    pts = vtk.vtkPoints()  
    ids = vtk.vtkIdList()  

    if tree.IntersectWithLine(p0, p1, pts, ids):  
        p0 = np.array(p0)  
        p1 = np.array(p1)  
        # 绕路径起点和终点创建圆形区域  
        buffer_p0 = Point(p0).buffer(radius)  
        buffer_p1 = Point(p1).buffer(radius)  
        
        for i in range(pts.GetNumberOfPoints()):  
            q = np.array(pts.GetPoint(i))  
            if buffer_p0.contains(Point(q)) or buffer_p1.contains(Point(q)):  
                return True  
    return False  

def make_tube(p0, p1):  
    pts = vtk.vtkPoints()   
    pts.InsertNextPoint(p0)   
    pts.InsertNextPoint(p1)  
    ln = vtk.vtkLine()   
    ln.GetPointIds().SetId(0, 0)   
    ln.GetPointIds().SetId(1, 1)  
    ca = vtk.vtkCellArray()   
    ca.InsertNextCell(ln)  
    pd = vtk.vtkPolyData()   
    pd.SetPoints(pts)   
    pd.SetLines(ca)  
    tf = vtk.vtkTubeFilter()   
    tf.SetInputData(pd)  
    tf.SetRadius(TUBE_DIAM / 2)   
    tf.SetNumberOfSides(20)   
    tf.CappingOn()   
    tf.Update()  
    return tf.GetOutput()  

def save_poly(pd, fname):  
    fn = os.path.join(BASE_DIR, fname)  
    wr = vtk.vtkPolyDataWriter()   
    wr.SetFileName(fn)  
    wr.SetInputData(pd)   
    wr.SetFileTypeToBinary()   
    wr.Write()  
    return fn  

# ───────────────────── ⑤ 读取模型 ─────────────────────  
print("加载模型 ...")  
entry_pd = load_poly(entry_path)  
hard_pd = load_poly(hard_path)  

# 提取入口点  
entry_pts = poly_points(entry_pd)[::PLANE_SKIP]  

targets = []  
for json_name, folder, color in ENDPOINTS:  
    js = os.path.join(BASE_DIR, json_name)  
    for pos in read_endpoints(js):  
        targets.append((folder, color, pos))  

print(f"入口点数: {len(entry_pts)}   目标点数: {len(targets)}")  

# ───────────────────── ⑥ 创建加速结构 ─────────────────────  
hard_tree = vtk.vtkOBBTree()  
hard_tree.SetDataSet(hard_pd)  
hard_tree.BuildLocator()  

# ───────────────────── ⑦ 分别计算每个终点的最短路径 ─────────────────────  
for folder, color, target in targets:  
    print(f"\n计算目标 {folder} 的路径...")  
    safe = []  # 存储无碰撞的安全路径  
    collisions = []  # 存储与障碍物相撞的路径  
    total = len(entry_pts)  
    bar = tqdm(total=total, ncols=80, unit='pair', bar_format='{l_bar}{bar}|')  

    for p0_raw in entry_pts:  
        bar.update(1)  
        
        v = np.array(target) - np.array(p0_raw)  
        L = np.linalg.norm(v)  
        if L < 1e-6: continue  
        
        dir_v = v / L  

        # 计算反向延长的起点  
        p0 = p0_raw - dir_v * EXTEND_LEN  
        Ltot = L + EXTEND_LEN  

        # 检查路径是否与障碍物相交  
        if hits_cylinder(p0, target, hard_tree, TUBE_DIAM / 2):  
            collisions.append((folder, color, tuple(p0), tuple(target), Ltot))  # 保存碰撞路径  
            continue  # 如果相交，跳过这个路径  

        # 记录有效路径  
        safe.append(dict(folder=folder, color=color,  
                         p0=tuple(p0), p1=tuple(target),  
                         length=Ltot))  

    bar.close()  
    print(f"\n有效路径总数: {len(safe)}")  
    print(f"相撞路径总数: {len(collisions)}")  

    # ───────────────────── ⑧ 保存有效路径 ─────────────────────  
    exports = []  

    # 输出最短路径  
    if safe:  
        min_length_path = min(safe, key=lambda p: p['length'])  
        print(f"到目标 {folder} 的最短路径: ")  
        print(f"    Length L={min_length_path['length']:.1f} mm")  
        poly = make_tube(min_length_path['p0'], min_length_path['p1'])  
        fn = save_poly(poly, f"{min_length_path['folder']}-ShortestPath.vtk")  
        exports.append((poly, min_length_path['color']))  
        print("    输出:", fn)  

# ───────────────────── ⑨ 可视化 ─────────────────────  
ren = vtk.vtkRenderer()  
ren.SetBackground(1, 1, 1)  

# 可视化所有有效路径  
for pd, col in exports:  
    mp = vtk.vtkPolyDataMapper()  
    mp.SetInputData(pd)  
    ac = vtk.vtkActor()  
    ac.SetMapper(mp)  
    ac.GetProperty().SetColor(*col)  
    ren.AddActor(ac)  

# 展示窗口  
win = vtk.vtkRenderWindow()  
win.AddRenderer(ren)  
win.SetSize(1100, 850)  
iren = vtk.vtkRenderWindowInteractor()  
iren.SetRenderWindow(win)  
ren.ResetCamera()  
win.Render()  
print("\n窗口已打开：显示所有最短路径")  
iren.Start()
