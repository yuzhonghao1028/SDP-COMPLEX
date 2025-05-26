import os  
import sys  
import json  
import time  
import numpy as np  
import vtk  
import tkinter as tk  
from tkinter import filedialog  
from tqdm import tqdm  
from shapely.geometry import Point  
import matplotlib.pyplot as pl

# ───────────────────── ① 输入用户文件夹 ─────────────────────  
def get_user_input():  
    root = tk.Tk()  
    root.withdraw()  

    base_dir = filedialog.askdirectory(title="请选择文件夹")  
    if not base_dir:  
        sys.exit("[Err] 文件夹未选择！")  

    return base_dir  

# ───────────────────── ② 固定参数 ─────────────────────  
EXTEND_LEN = 20.0  
IGNORE_R = 3.0  
TUBE_DIAM = 3.0  # 圆柱体直径为3mm
  
# ───────────────────── ③ 固定目录 ─────────────────────  
BASE_DIR = get_user_input()  
entry_path  = os.path.join(BASE_DIR, "entry.vtk")  
hard_path   = os.path.join(BASE_DIR, "no surgery area.vtk")  

ENDPOINTS = [  
    ("anterior.mrk.json",  "anterior"),   
    ("middle.mrk.json",    "middle"),  
    ("posterior.mrk.json", "superior"),   
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

def save_poly(pd, fname):  
    fn = os.path.join(BASE_DIR, fname)  
    wr = vtk.vtkPolyDataWriter()   
    wr.SetFileName(fn)  
    wr.SetInputData(pd)   
    wr.SetFileTypeToBinary()   
    wr.Write()  
    return fn  

def point_to_line_distance(point, line_start, line_end):  
    """计算点到线段的距离"""  
    point = np.array(point)  
    line_start = np.array(line_start)  
    line_end = np.array(line_end)  

    # 计算线段向量和点到线段起点的向量  
    line_vec = line_end - line_start  
    point_vec = point - line_start  

    # 计算线段长度的平方  
    line_len_squared = np.dot(line_vec, line_vec)  
    if line_len_squared == 0:  
        return np.linalg.norm(point - line_start)  # 线段长度为零，返回起点到点的距离  

    # 计算投影长度  
    projection_length = np.dot(point_vec, line_vec) / line_len_squared  

    # 限制投影长度在[0, 1]范围内  
    if projection_length < 0:  
        closest_point = line_start  
    elif projection_length > 1:  
        closest_point = line_end  
    else:  
        closest_point = line_start + projection_length * line_vec  

    return np.linalg.norm(point - closest_point)  # 返回点到最近点的距离  

def create_line(start, end):  
    """创建连接起点和终点的线段"""  
    line_points = vtk.vtkPoints()  
    line_points.InsertNextPoint(start)  
    line_points.InsertNextPoint(end)  

    line = vtk.vtkCellArray()  
    line.InsertNextCell(2)  # Line has 2 points  
    line.InsertCellPoint(0)  
    line.InsertCellPoint(1)  

    line_polydata = vtk.vtkPolyData()  
    line_polydata.SetPoints(line_points)  
    line_polydata.SetLines(line)  

    return line_polydata


# ───────────────────── ⑤ 读取模型 ─────────────────────  
print("加载模型 ...")  
entry_pd = load_poly(entry_path)  
hard_pd = load_poly(hard_path)  

# 提取入口点
entry_pts = poly_points(entry_pd)[::2]  

targets = []  
for json_name, folder in ENDPOINTS:  # 这里只解包两个元素  
    js = os.path.join(BASE_DIR, json_name)  
    for pos in read_endpoints(js):  
        targets.append((folder, pos))  # 直接使用folder和pos

print(f"入口点数: {len(entry_pts)}   目标点数: {len(targets)}")

# ───────────────────── ⑥ 计算有效路径 ─────────────────────  
final_paths = []  # 存储每个目标的有效路径信息  

for folder, target in targets:  # 这里只解包 folder 和 target  
    print(f"\n计算目标 {folder} 的有效路径...")  
    safe_paths = []  # 存储有效路径的列表  

    # 使用 tqdm 包装 entry_pts 以显示进度条  
    for p0_raw in tqdm(entry_pts, desc=f"计算 {folder} 的路径", unit="起点"):  
        collision_detected = False  
        
        model_points = poly_points(hard_pd)  

        for model_point in model_points:  
            distance = point_to_line_distance(model_point, p0_raw, target)  
            if distance <= IGNORE_R:  
                collision_detected = True  
                break  

        if not collision_detected:  
            path_length = np.linalg.norm(np.array(target) - np.array(p0_raw))  # 计算路径长度  
            safe_paths.append((tuple(p0_raw), target, path_length))  # 确保加入路径长度  

if safe_paths:  
    min_length_path = min(safe_paths, key=lambda p: p[2])  # 找到长度最短的路径  
    # 只保存起点和终点，不再包括颜色  
    final_paths.append((min_length_path[0], min_length_path[1]))  
    print(f"有效路径到目标 {folder}：从 {min_length_path[0]} 到 {min_length_path[1]}，长度 {min_length_path[2]:.2f} mm")  
else:  
    print(f"没有有效路径到达目标 {folder}.")
# ───────────────────── ⑦ 创建并保存圆柱体与线段 ─────────────────────  
if final_paths:  
    for index, (folder, _) in enumerate(ENDPOINTS):  
        if index < len(final_paths):  
            start, end = final_paths[index]  
            target_name = f"{folder} path"  
            
            # 创建连接线段  
            line = create_line(start, end)  

            # 保存线段   
            fn_line = save_poly(line, f"{target_name}_line.vtk")  
            print(f"输出: {fn_line}")  
        else:  
            print(f"警告: 对于 {folder}，无效路径。没有找到对应的路径。")  
else:  
    print("没有有效路径可供绘制！")  
