import numpy as np
from scipy.spatial import Voronoi
from scipy.spatial import cKDTree
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator

def parse_nodes_and_elements_fast(file_path):
    # 使用字典存储：ID 作为键，数据作为值
    node_dict = {}    # {ID: [x, y, z]}
    element_dict = {} # {type: {ID: [n1, n2, ...]}}
    
    current_section = None
    current_el_type = None

    with open(file_path, 'r') as f:
        for line in f:
            line_raw = line.strip()
            if not line_raw or line_raw.startswith('**'):
                continue
            
            line_compact = line_raw.upper().replace(" ", "")
            
            if line_compact.startswith('*'):
                if line_compact == '*NODE':
                    current_section = 'NODE'
                    continue
                elif line_compact.startswith('*ELEMENT,TYPE='):
                    current_section = 'ELEMENT'
                    current_el_type = line_compact.split('=')[1]
                    if current_el_type not in element_dict:
                        element_dict[current_el_type] = {}
                    continue
                else:
                    current_section = None
                    continue

            # --- 解析数据行 ---
            parts = line_raw.split(',')
            if current_section == 'NODE':
                try:
                    # parts[0] 是节点号，parts[1:] 是坐标
                    node_id = int(float(parts[0]))
                    coords = [float(x) for x in parts[1:]]
                    node_dict[node_id] = coords
                except (ValueError, IndexError):
                    continue
                    
            elif current_section == 'ELEMENT':
                try:
                    # parts[0] 是单元号，parts[1:] 是连接的节点号
                    el_id = int(float(parts[0]))
                    connectivity = [int(float(x)) for x in parts[1:] if x.strip()]
                    element_dict[current_el_type][el_id] = connectivity
                except (ValueError, IndexError):
                    continue

    return node_dict, element_dict
    nodes = []
    elements = {} 
    
    current_section = None
    current_el_type = None

    print(f"正在解析文件: {file_path} ...")
    
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line_raw = line.strip()
            # 跳过空行和注释
            if not line_raw or line_raw.startswith('**'):
                continue
            
            # 将行转为大写并去除所有空格，用于严格匹配
            line_compact = line_raw.upper().replace(" ", "")
            
            # --- 严格匹配关键词行 ---
            if line_compact.startswith('*'):
                if line_compact == '*NODE':
                    current_section = 'NODE'
                    continue
                
                elif line_compact == '*ELEMENT,TYPE=C3D8':
                    current_section = 'ELEMENT'
                    current_el_type = 'C3D8'
                    if current_el_type not in elements:
                        elements[current_el_type] = []
                    continue
                
                elif line_compact == '*ELEMENT,TYPE=C3D6':
                    current_section = 'ELEMENT'
                    current_el_type = 'C3D6'
                    if current_el_type not in elements:
                        elements[current_el_type] = []
                    continue
                
                else:
                    # 碰到其他任何关键词（如 *OUTPUT, *STEP 等）则停止记录数据
                    current_section = None
                    continue

            # --- 解析数据行 ---
            if current_section == 'NODE':
                try:
                    # 存储格式: [节点号, x, y, z]
                    nodes.append([float(x) for x in line_raw.split(',')])
                except ValueError:
                    continue
                    
            elif current_section == 'ELEMENT':
                try:
                    # 存储格式: [单元号, 节点1, 节点2, ...]
                    # 注意：Abaqus 单元节点可能通过逗号换行，这里默认处理单行完整数据
                    elements[current_el_type].append([int(float(x)) for x in line_raw.split(',')])
                except ValueError:
                    continue

    # 转换为 Numpy 数组
    node_array = np.array(nodes)
    element_arrays = {k: np.array(v) for k, v in elements.items()}
    
    return node_array, element_arrays

# !!!检查自适应网格Adaptive Mesh的节点集和单元集名称
def parse_sets_from_inp(file_path):

    """
    专门解析指定名称的 Nset 和 Elset。
    逻辑：先定位 User Constraint 关键字，获取 Set 名称，再解析对应数据。
    """
    target_set_name = None
    node_set_data = set()
    elem_set_data = set()
    
    # 第一次读取：快速定位 *Adaptive Mesh Constraint 下的名称
    with open(file_path, 'r') as f:
        lines = f.readlines() # 如果文件极大（GB级），建议改用迭代器，这里假设内存充足
        for i, line in enumerate(lines):
            if '*ADAPTIVE MESH CONSTRAINT, USER' in line.upper():
                # 获取下一行作为目标名称
                target_set_name = lines[i+1].strip()
                break
    
    if not target_set_name:
        print("未找到 Adaptive Mesh Constraint 指定的 Set 名称。")
        return None

    print(f"识别到目标 Set 名称: {target_set_name}")

    # 第二次解析：提取 Nset 和 Elset 数据
    current_section = None
    is_generate = False
    
    # 构造识别头
    nset_header = f"*NSET,NSET={target_set_name.upper()}"
    elset_header = f"*ELSET,ELSET={target_set_name.upper()}"
    ####elset_header = f"*ELSET,ELSET={'_PickedSet134'.upper()}" #非abaqus划分网格的inp文件

    for line in lines:
        line_raw = line.strip()
        if not line_raw or line_raw.startswith('**'):
            continue
            
        line_upper_compact = line_raw.upper().replace(" ", "")

        # 检查是否进入目标区块
        if line_upper_compact.startswith('*'):
            is_generate = 'GENERATE' in line_upper_compact
            if nset_header in line_upper_compact:
                current_section = 'NSET'
                continue
            elif elset_header in line_upper_compact:
                current_section = 'ELSET'
                continue
            else:
                current_section = None # 碰到其他关键词，结束当前区块
                continue

        # 解析数据行
        if current_section in ['NSET', 'ELSET']:
            try:
                # 提取数字
                values = [int(float(x)) for x in line_raw.strip(',').split(',')]
                
                target_container = node_set_data if current_section == 'NSET' else elem_set_data
                
                if is_generate:
                    # 如果有 generate 关键字：起始, 终止, 步长
                    if len(values) >= 2:
                        start = values[0]
                        end = values[1]
                        step = values[2] if len(values) > 2 else 1
                        # 使用 range 快速生成并更新集合
                        target_container.update(range(start, end + 1, step))
                else:
                    # 普通列表格式：直接添加所有 ID
                    target_container.update(values)
            except ValueError:
                continue

    # 返回结果转为 numpy 数组方便后续计算，或保持 set 方便匹配
    return {
        'name': target_set_name,
        'nodes': np.sort(list(node_set_data)),
        'elements': np.sort(list(elem_set_data))
    }

    """
    在指定的节点集合中，筛选出 Z 坐标在指定范围内的节点。
    
    参数:
    - nodes_dict: 全量节点字典 {ID: [x, y, z]}
    - target_node_ids: 待筛选的节点 ID 集合 (Set)
    - z_target: 目标 Z 坐标
    - tolerance: 容差范围
    
    返回:
    - filtered_nodes: 筛选后的节点字典 {ID: [x, y, z]}
    """
    # 计算上下界，避免在循环中重复计算
    lower_bound = z_target - tolerance
    upper_bound = z_target + tolerance
    
    # 使用字典推导式高效筛选
    # 逻辑：只遍历 target_node_ids（缩小搜索范围），利用 nodes_dict 的 O(1) 查找获取坐标
    filtered_nodes = {
        node_id: nodes_dict[node_id] 
        for node_id in target_node_ids 
        if node_id in nodes_dict and (lower_bound <= nodes_dict[node_id][2] <= upper_bound)
    }
    
    return filtered_nodes

def filter_element_nodes_by_z(nodes_dict, element_dict, set_results, z_target=0.0, tolerance=10):
    """
    根据Z坐标筛选单元关联的节点。
    nodes_dict: {node_id: [x, y, z]}
    element_dict: {type: {el_id: [n1, n2, ...]}}
    set_results: {'elements': [id1, id2, ...], 'nodes': [...]}
    """
    # 1. 预计算：在整个 node_dict 中，找出所有满足 Z 坐标条件的节点 ID，存入 set
    # 这一步是为了将百万次的坐标判断简化为一次 O(1) 的集合查找
    valid_z_node_ids = set()
    for node_id, coords in nodes_dict.items():
        # coords[2] 是 Z 坐标
        if abs(coords[2] - z_target) <= tolerance:
            valid_z_node_ids.add(node_id)
    
    # 2. 获取目标单元 ID 集合（转为 set 提高匹配速度）
    target_el_ids = set(set_results['elements'])
    
    # 3. 结果存储
    filtered_elements = {}

    # 4. 遍历所有单元类型
    for el_type, el_data in element_dict.items():
        filtered_elements[el_type] = {}
        
        for el_id, connectivity in el_data.items():
            # 只处理在 set_results 指定集合中的单元
            ####if el_id==3106009:
            ####    print("debug")
            if el_id in target_el_ids:
                # 核心过滤逻辑：保留在 valid_z_node_ids 中的节点
                # 使用列表推导式，这是 Python 处理此类逻辑最快的方式
                new_conn = [n_id for n_id in connectivity if n_id in valid_z_node_ids]
                
                # 如果过滤后该单元还有节点剩余，则保留该单元
                if new_conn:
                    filtered_elements[el_type][el_id] = new_conn
            else:
                # 如果单元不在指定的 set 中，是否保留原始数据？
                # 这里假设你只想获取过滤后的集合结果，如果需要保留全部，请取消下行注释
                # filtered_elements[el_type][el_id] = connectivity
                pass

    return filtered_elements

def build_node_adjacency(filtered_el_dict, target_node_set):
    """
    构建节点邻接表：找出指定节点集中的每个点所连接的其他点。
    filtered_el_dict: {type: {el_id: [n1, n2, ...]}}
    target_node_set: set 或 list，包含需要查询的目标节点ID
    """
    # 将目标节点转为 set，实现 O(1) 查询
    nodes_of_interest = set(target_node_set)
    
    # 初始化邻接表：{节点ID: {邻居1, 邻居2, ...}}
    # 使用 set 存储邻居，自动去重
    adjacency_map = {node_id: set() for node_id in nodes_of_interest}

    print("开始构建节点邻接关系...")

    # 遍历所有单元类型和数据
    for el_type, el_map in filtered_el_dict.items():
        for connectivity in el_map.values():
            # 过滤出当前单元中属于目标节点集的节点
            # 只有这些节点才需要记录它们的邻居
            current_nodes = [n for n in connectivity if n in nodes_of_interest]
            
            if len(current_nodes) < 2:
                continue  # 如果单元内目标节点少于2个，无法构成连接
            
            # 单元内节点两两互联
            # 这种写法比嵌套循环更快：每个节点连接单元内除自己外的所有节点
            unit_nodes_set = set(connectivity) 
            for node_id in current_nodes:
                # 将当前单元的所有节点加入该节点的邻居集
                adjacency_map[node_id].update(unit_nodes_set)

    # 最后处理：移除节点自身，并根据需要转换为列表
    final_adjacency = {}
    for node_id, neighbors in adjacency_map.items():
        if node_id in neighbors:
            neighbors.remove(node_id) # 移除自己连接自己
        
        # 如果需要有序或者方便后续处理，可以转为 np.array 或 list
        final_adjacency[node_id] = np.array(list(neighbors), dtype=int)

    return final_adjacency

def sort_node_neighbors_custom(node_neighbors, nodes_dict):
    """
    根据特定规则对邻接节点进行排序。
    node_neighbors: {node_id: np.array([neighbor_ids])}
    nodes_dict: {node_id: [x, y, z]}
    """
    sorted_neighbors = {}

    def get_dist_sq(id1, id2):
        # 使用平方距离比较，避免计算开方(sqrt)，速度更快
        p1 = nodes_dict[id1]
        p2 = nodes_dict[id2]
        return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2

    for main_id, neighbors in node_neighbors.items():
        n_count = len(neighbors)
        if n_count == 0:
            sorted_neighbors[main_id] = neighbors
            continue
            
        main_coords = nodes_dict[main_id]
        
        # --- 规则 1: 数量等于 3 ---
        if n_count == 3:
            # 计算到主节点的距离
            dists = [(nid, get_dist_sq(main_id, nid)) for nid in neighbors]
            # 距离最远的排最后
            dists.sort(key=lambda x: x[1])
            farthest_node = dists[-1][0]
            remains = [dists[0][0], dists[1][0]]
            # 剩余节点按 |x - x_main| 排序
            remains.sort(key=lambda nid: abs(nodes_dict[nid][0] - main_coords[0]), reverse=True)
            sorted_neighbors[main_id] = np.array([remains[0], remains[1], farthest_node])

        # --- 规则 2: 数量等于 5 ---
        elif n_count == 5:
            dists = [(nid, get_dist_sq(main_id, nid)) for nid in neighbors]
            dists.sort(key=lambda x: x[1])
            # 距离最远两个
            farthest_two = [dists[-2][0], dists[-1][0]]
            remains = [dists[0][0], dists[1][0], dists[2][0]]
            # 剩余按 (main_x - node_x) 排序：最大、最小、中间
            remains.sort(key=lambda nid: main_coords[0] - nodes_dict[nid][0], reverse=True)
            sorted_neighbors[main_id] = np.array([remains[0], remains[2], remains[1], farthest_two[0], farthest_two[1]])

        # --- 规则 3: 数量大于 5 ---
        elif n_count > 5:
            # 提取所有邻居坐标
            neighbor_coords = np.array([nodes_dict[nid] for nid in neighbors])
            xs = neighbor_coords[:, 0]
            ys = neighbor_coords[:, 1]
            
            xmin, xmax = np.min(xs), np.max(xs)
            ymin, ymax = np.min(ys), np.max(ys)
            y_mid = (ymax - ymin) / 2 + ymin
            x_mid = (xmax - xmin) / 2 + xmin
            
            # 4个虚拟点
            v_points = [
                np.array([xmin, y_mid]),  # V1
                np.array([xmax, y_mid]),  # V2
                np.array([x_mid, ymax]),  # V3
                np.array([x_mid, ymin])   # V4
            ]
            
            # 找出最靠近4个虚拟点的节点
            assigned_ids = []
            pool = list(neighbors)
            
            for vp in v_points:
                # 计算pool中所有点到当前虚拟点的2D距离
                best_node = min(pool, key=lambda nid: (nodes_dict[nid][0]-vp[0])**2 + (nodes_dict[nid][1]-vp[1])**2)
                assigned_ids.append(best_node)
                pool.remove(best_node)
            
            # 剩余节点按到主节点的直线距离排序
            pool.sort(key=lambda nid: get_dist_sq(main_id, nid))
            
            sorted_neighbors[main_id] = np.array(assigned_ids + pool)

        else:
            # 其他数量（1, 2, 4）默认按距离排序或保持原样
            sorted_neighbors[main_id] = neighbors

    return sorted_neighbors


    """
    计算主节点与其邻接节点围成的特定面积
    sorted_neighbors: {main_id: np.array([ordered_neighbors])}
    nodes_dict: {node_id: [x, y, z]}
    """
    node_areas = {}

    def get_poly_area_2d(coords):
        """使用鞋带公式计算2D多边形面积: 0.5 * |Σ(xi*yi+1 - xi+1*yi)|"""
        x = coords[:, 0]
        y = coords[:, 1]
        return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

    print("开始计算几何面积...")

    for main_id, neighbors in sorted_neighbors.items():
        n_count = len(neighbors)
        if n_count < 3:
            node_areas[main_id] = 0.0
            continue

        # 获取主节点坐标 (取 X, Y 投影面进行面积计算)
        main_pt = np.array(nodes_dict[main_id][:2])
        
        # --- 规则 1: 3个邻接点 (主节点+3个邻居=4点四边形) ---
        if n_count == 3:
            # 组合 4 个点的坐标
            pts = np.array([main_pt] + [nodes_dict[nid][:2] for nid in neighbors])
            # 注意：此处假设点集已按顺/逆时针排列，若非凸四边形需重新排序
            node_areas[main_id] = get_poly_area_2d(pts)

        # --- 规则 2: 5个邻接点 (前3个邻居形成的三角形) ---
        elif n_count == 5:
            pts = np.array([nodes_dict[nid][:2] for nid in neighbors[:3]])
            node_areas[main_id] = get_poly_area_2d(pts)

        # --- 规则 3: >5个邻接点 (主节点与前4个邻居的特征面积) ---
        elif n_count > 5:
            # 根据前述排序，前4个点分别是靠近 xmin, xmax, ymax, ymin 的点
            # 这四个点构成的四边形覆盖了主节点的核心受载区域（Voronoi 面积的快速估算）
            pts = np.array([nodes_dict[nid][:2] for nid in neighbors[:4]])
            node_areas[main_id] = get_poly_area_2d(pts)
            
        else:
            node_areas[main_id] = 0.0

    return node_areas

def calculate_neighbor_areas(sorted_neighbors, nodes_dict):
    """
    计算主节点与其邻接节点围成的特定面积
    针对非排序点集进行稳健面积计算
    """
    node_areas = {}

    def get_tri_area_2d(p1, p2, p3):
        """计算三个点组成的三角形面积"""
        # 向量叉乘公式: 0.5 * |x1(y2-y3) + x2(y3-y1) + x3(y1-y2)|
        return 0.5 * abs(p1[0]*(p2[1]-p3[1]) + p2[0]*(p3[1]-p1[1]) + p3[0]*(p1[1]-p2[1]))

    print("开始计算几何面积...")

    for main_id, neighbors in sorted_neighbors.items():
        n_count = len(neighbors)
        if n_count < 3:
            node_areas[main_id] = 0.0
            continue

        # 获取主节点坐标 (X, Y)
        m = np.array(nodes_dict[main_id][:2])
        
        # --- 规则 1: 3个邻接点 (主节点+3个邻居) ---
        if n_count == 3:
            # 设邻居为 n0, n1, n2。主节点为 m。
            # 为了处理非顺时针点集，最稳健的方法是计算：
            # 以主节点为中心，与邻居两两构成的三角形面积之和
            n0 = np.array(nodes_dict[neighbors[0]][:2])
            n1 = np.array(nodes_dict[neighbors[1]][:2])
            n2 = np.array(nodes_dict[neighbors[2]][:2])
            
            # 这里根据物理含义，通常主节点在中间，面积 = Area(m,n0,n1) + Area(m,n1,n2) + Area(m,n2,n0)
            # 这种方法不依赖邻居的全局排序，只依赖它们绕主节点的局部位置
            area = get_tri_area_2d(m, n1, n2) + get_tri_area_2d(m, n2, n0)
            node_areas[main_id] = area

        # --- 规则 2: 5个邻接点 (前3个邻居形成的三角形) ---
        elif n_count == 5:
            # 前3个邻居按照之前的规则已选出，直接计算这三点三角形
            n0 = np.array(nodes_dict[neighbors[0]][:2])
            n1 = np.array(nodes_dict[neighbors[1]][:2])
            n2 = np.array(nodes_dict[neighbors[2]][:2])
            node_areas[main_id] = get_tri_area_2d(n0, n1, n2)

        # --- 规则 3: >5个邻接点 (主节点与前4个邻居) ---
        elif n_count > 5:
            try:
                # 仅取主节点 + 前4个邻居（共5个点）
                # 这4个点已按之前的虚拟点规则排序，分布在四个象限/方位
                pts = [m] + [nodes_dict[neighbors[i]][:2] for i in range(4)]
                pts = np.array(pts)
                
                # 构建 5 个点的微型 Voronoi
                vor = Voronoi(pts)
                
                # 获取主节点（索引0）的区域
                region_idx = vor.point_region[0]
                vertices_indices = vor.regions[region_idx]
                
                # 判定：如果有无限区域(-1)或数据不全，使用三角形保底
                if -1 in vertices_indices or not vertices_indices:
                    p = pts[1:] # 4个邻居
                    area = get_tri_area_2d(m, p[0], p[1]) + get_tri_area_2d(m, p[1], p[2]) + \
                           get_tri_area_2d(m, p[2], p[3]) + get_tri_area_2d(m, p[3], p[0])
                else:
                    # 鞋带公式计算有限多边形面积
                    v = vor.vertices[vertices_indices]
                    area = 0.5 * np.abs(np.dot(v[:,0], np.roll(v[:,1], 1)) - np.dot(v[:,1], np.roll(v[:,0], 1)))
                
                node_areas[main_id] = area
            except:
                node_areas[main_id] = 0.0
            
        else:
            node_areas[main_id] = 0.0

    return node_areas

def identify_node_attributes(nodes_dict, target_node_set, tolerance=1e-5):
    """
    判别节点集中每个节点的边界属性。
    1: ymin, 3: ymax, 2: xmax, 4: xmin, 0: 内部点
    """
    # 1. 提取目标节点坐标到 Numpy 数组以实现极速计算
    target_node_ids = np.array(list(target_node_set))
    # 提取对应的 (x, y) 坐标
    coords = np.array([nodes_dict[nid][:2] for nid in target_node_ids])
    
    xs = coords[:, 0]
    ys = coords[:, 1]

    # 2. 计算全局极值
    xmin, xmax = np.min(xs), np.max(xs)
    ymin, ymax = np.min(ys), np.max(ys)
    
    print(f"坐标范围: X[{xmin}, {xmax}], Y[{ymin}, {ymax}]")

    # 3. 按照顺序进行属性判别（使用 Numpy 掩码/Mask 提高效率）
    # 初始化所有属性为 0
    attributes = np.zeros(len(target_node_ids), dtype=int)
    
    # 优先级 1: ymin (属性 1)
    mask1 = np.abs(ys - ymin) <= tolerance
    attributes[mask1] = 1
    
    # 优先级 2: ymax (属性 3，且之前未被定义)
    mask3 = (np.abs(ys - ymax) <= tolerance) & (attributes == 0)
    attributes[mask3] = 3
    
    # 优先级 3: xmax (属性 2，且之前未被定义)
    mask2 = (np.abs(xs - xmax) <= tolerance) & (attributes == 0)
    attributes[mask2] = 2
    
    # 优先级 4: xmin (属性 4，且之前未被定义)
    mask4 = (np.abs(xs - xmin) <= tolerance) & (attributes == 0)
    attributes[mask4] = 4

    # 4. 将结果映射回字典 {node_id: attribute}
    node_attr_map = dict(zip(target_node_ids, attributes))
    
    return node_attr_map

def interpolate_precipitation(nodes_dict, target_node_set, prec_file_path):
    """
    基于 precipitation.txt 的 (x, y, pre) 数据，为目标节点插值 pre 值。
    使用 cKDTree 实现百万级数据的极速最近邻插值。
    """
    print(f"开始读取降水数据: {prec_file_path} ...")
    
    # 1. 读取降水数据 (txt 文件通常为：x, y, pre)
    # 使用 np.loadtxt 并指定建议的 fast 引擎（如果可用）
    try:
        prec_data = np.loadtxt(prec_file_path, delimiter=',') # 如果是空格分隔，去掉 delimiter
    except ValueError:
        prec_data = np.loadtxt(prec_file_path) # 尝试默认空格/制表符分隔
        
    prec_coords = prec_data[:, :2] # 降水观测点的 x, y
    prec_values = prec_data[:, 2]  # 对应的降水值 pre

    # 2. 构建 k-d Tree (空间索引)
    print("构建空间索引树...")
    tree = cKDTree(prec_coords)

    # 3. 准备目标节点的坐标
    target_node_ids = np.array(list(target_node_set))
    # 提取目标节点的 (x, y)
    target_coords = np.array([nodes_dict[nid][:2] for nid in target_node_ids])

    # 4. 执行查询：寻找每个节点最近的降水观测点
    # k=1 表示最近邻插值；如果需要更平滑的线性插值，可设 k=4 并计算加权平均
    print(f"正在为 {len(target_node_ids)} 个节点进行空间插值...")
    distances, indices = tree.query(target_coords, k=1)

    # 5. 获取插值结果
    # 根据索引直接从 prec_values 中提取对应的 pre 值
    interpolated_pre = prec_values[indices]

    # 6. 将结果映射回字典 {node_id: pre_value}
    node_pre_map = dict(zip(target_node_ids, interpolated_pre))

    print("插值完成。")
    return node_pre_map

def reindex_and_format_neighbors(target_node_set, sorted_neighbors):
    """
    1. 对目标节点集重新编号 (1, 2, 3...)
    2. 更新 sorted_neighbors 中的邻接节点 ID
    3. 统一格式化为 8 个数值（截断或补0）
    """
    # --- 第一步：建立 ID 映射表 ---
    # 将 set 转换为有序列表（从小到大）
    sorted_original_ids = sorted(list(target_node_set))
    
    # 建立旧 ID 到新 ID 的字典映射 {old_id: new_id}
    # 使用 dict 确保百万级查询速度为 O(1)
    id_mapping = {old_id: i + 1 for i, old_id in enumerate(sorted_original_ids)}
    
    # 为了处理邻接节点中可能存在的不在 target_node_set 里的节点
    # 我们给 mapping 增加一个默认行为：如果没找到，映射为 0
    # 或者根据你的需求，不在 set 里的点也计为 0
    
    final_reindexed_neighbors = {}

    print(f"开始重编号与格式化，涉及主节点数: {len(sorted_neighbors)} ...")

    # --- 第二步：遍历主节点并转换邻居 ID ---
    for main_id, neighbors in sorted_neighbors.items():
        # 1. 将邻居旧 ID 转换为新 ID（如果邻居不在 mapping 中，则转换为 0）
        new_neighbor_ids = [id_mapping.get(nid, 0) for nid in neighbors]
        
        # 2. 截断或补齐至 8 个数值
        if len(new_neighbor_ids) >= 8:
            # 只取前 8 个
            formatted_list = new_neighbor_ids[:8]
        else:
            # 在后面补 0 直至 8 个
            formatted_list = new_neighbor_ids + [0] * (8 - len(new_neighbor_ids))
        
        # 3. 存入结果字典
        final_reindexed_neighbors[main_id] = formatted_list

    print("重编号与格式化完成。")
    return final_reindexed_neighbors, id_mapping

def output_model_files(final_map, node_areas, node_attributes, nodes_dict, node_pre_results):
    """
    输出 ada_ele.txt, mo_pre.txt 和 model_param.txt。
    支持百万级数据，使用流式写入以节省内存。
    """
    
    # 按照主节点编号从小到大进行排序输出，方便后续匹配
    sorted_main_ids = sorted(final_map.keys())
    
    # --- 1. 输出 ada_ele.txt ---
    print("正在生成 ada_ele.txt ...")
    with open('ada_ele.txt', 'w', encoding='utf-8') as f:
        for mid in sorted_main_ids:
            # 邻接节点列表 (5-12列)
            neighbor_list = final_map[mid]
            # 计算不为0的邻接节点总数
            count_nonzero = sum(1 for x in neighbor_list if x != 0)
            # 获取面积和属性
            area = node_areas.get(mid, 0.0)
            attr = node_attributes.get(mid, 0)
            
            # 组合成一行：ID, Area, Attr, Count, N1, N2, N3, N4, N5, N6, N7, N8
            neighbors_str = " ".join(map(str, neighbor_list))
            f.write(f"{mid} {round(area)} {attr} {count_nonzero} {neighbors_str}\n")

    # --- 2. 输出 mo_pre.txt ---
    print("正在生成 mo_pre.txt ...")
    with open('mo_pre.txt', 'w', encoding='utf-8') as f:
        for mid in sorted_main_ids:
            # 获取坐标
            coords = nodes_dict.get(mid, [0.0, 0.0, 0.0])
            # 获取降雨插值，如果为空则为 0
            pre = node_pre_results.get(mid, 0.0)
            
            # 组合成一行：ID, X, Y, Z, Pre
            f.write(f"{mid} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f} {pre:.6f}\n")

    # --- 3. 输出 model_param.txt ---
    print("正在生成 model_param.txt ...")
    if sorted_main_ids:
        min_id = sorted_main_ids[0]
        max_id = sorted_main_ids[-1]
        total_count = len(sorted_main_ids)
        with open('model_param.txt', 'w', encoding='utf-8') as f:
            f.write(f"{min_id} {max_id} {total_count} {12}\n")
    
    print("所有文件输出完成。")

def update_nodes_z_with_nc(file_path_nc, nodes_dict, target_node_set, inp_old, inp_new):
    """
    1. 读取NC地形数据
    2. 建立经纬度与模型X,Y的映射
    3. 插值获取新Z坐标
    4. 快速生成替换Z后的新INP文件
    """
    # --- 1. 读取NC数据 ---
    print(f"正在读取 NC 文件: {file_path_nc}")
    nc_data = Dataset(file_path_nc)
    # 自动识别变量名
    all_vars = nc_data.variables.keys()
    
    # 识别经度
    lon_name = [v for v in all_vars if v.lower() in ['lon', 'longitude', 'x']][0]
    # 识别纬度
    lat_name = [v for v in all_vars if v.lower() in ['lat', 'latitude', 'y']][0]
    # 识别地形数据 (排除经纬度后的主变量)
    elev_name = [v for v in all_vars if v not in [lon_name, lat_name] and len(nc_data.variables[v].shape) >= 2][0]

    print(f"检测到变量名: 经度={lon_name}, 纬度={lat_name}, 地形={elev_name}")

    lon_arr = nc_data.variables[lon_name][:]
    lat_arr = nc_data.variables[lat_name][:]
    elev_arr = nc_data.variables[elev_name][:]
    
    # --- 2. 建立坐标映射关系 ---
    # 获取目标节点集的模型坐标极值
    target_node_ids = np.array(list(target_node_set))
    target_coords = np.array([nodes_dict[nid][:2] for nid in target_node_ids])
    x_min, x_max = np.min(target_coords[:, 0]), np.max(target_coords[:, 1])
    y_min, y_max = np.min(target_coords[:, 1]), np.max(target_coords[:, 1])
    
    # NC数据的经纬度范围
    lon_min, lon_max = 100.0, 107.0
    lat_min, lat_max = 28.0, 34.0

    def map_coords_to_geo(x, y):
        """将模型x,y映射到nc的经纬度"""
        lon = (x - x_min) / (x_max - x_min) * (lon_max - lon_min) + lon_min
        lat = (y - y_min) / (y_max - y_min) * (lat_max - lat_min) + lat_min
        return lon, lat

    # --- 3. 执行插值 ---
    print("正在构建网格插值器...")
    # RegularGridInterpolator 要求输入维度顺序与数据一致
    # 假设 elev_arr 索引是 [lat, lon]
    interp_func = RegularGridInterpolator((lat_arr, lon_arr), elev_arr, 
                                            bounds_error=False, fill_value=None)

    new_z_map = {}
    print(f"正在为 {len(target_node_ids)} 个节点计算新高度...")
    
    # 批量映射并插值
    target_lons, target_lats = map_coords_to_geo(target_coords[:, 0], target_coords[:, 1])
    # 准备插值点数据 [lat, lon]
    interp_pts = np.vstack([target_lats, target_lons]).T
    new_zs = interp_func(interp_pts)
    
    # 存入字典映射 {ID: new_z}
    new_z_map = dict(zip(target_node_ids, new_zs))

    # --- 4. 生成新 INP 文件 (快速流式读写) ---
    print(f"正在生成新文件: {inp_new}")
    current_section = None
    
    with open(inp_old, 'r',) as f_in, \
         open(inp_new, 'w') as f_out:
        
        for line in f_in:
            line_raw = line.strip()
            line_upper = line_raw.upper().replace(" ", "")
            
            if line_upper.startswith('*'):
                if line_upper == '*NODE':
                    current_section = 'NODE'
                else:
                    current_section = None
                f_out.write(line)
                continue
            
            if current_section == 'NODE':
                try:
                    parts = line_raw.split(',')
                    node_id = int(float(parts[0]))
                    if node_id in new_z_map:
                        # 替换 Z 坐标 (parts[3])
                        new_z = new_z_map[node_id]
                        # 保持原始格式：ID, X, Y, NewZ
                        new_line = f"{parts[0]}, {parts[1]}, {parts[2]}, {new_z: .6f}\n"
                        f_out.write(new_line)
                    else:
                        f_out.write(line)
                except:
                    f_out.write(line)
            else:
                f_out.write(line)

    print("地形同步完成。")

# !!!检查模型边界名称
def process_velocity_boundary(inp_path, velocity_file, nodes_dict):
    """
    1. 解析 Boundary 下的 Set 名称
    2. 提取 Nset 中的节点 ID
    3. 插值获取速度 (vx, vy)
    4. 导出 modelve.txt
    """
    boundary_set_name = None
    boundary_node_ids = set()

    # --- 1. 第一次扫描：寻找 *Boundary 下的 Set 名称 ---
    print("正在定位 Boundary User Set...")
    with open(inp_path, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            line_upper = line.upper().replace(" ", "")
            if '*BOUNDARY,USER' in line_upper or '*BOUNDARY,OP=NEW,USER' in line_upper:
            ####if '*BOUNDARY,OP=NEW,USER' in line_upper:    
                # 读取下一行获取 Set 名称 (例如 Set-31, 1, 1, 1. -> 提取 Set-31)
                next_line = lines[i+1].strip()
                if next_line and not next_line.startswith('*'):
                    boundary_set_name = next_line.split(',')[0].strip()
                    break
    
    if not boundary_set_name:
        print("错误：未找到 *Boundary, op=NEW, user 对应的 Set。")
        return

    print(f"识别到边界 Set 名称: {boundary_set_name}")

    # --- 2. 第二次扫描：提取 Nset 中的节点 ID ---
    current_section = None
    target_header = f"*NSET,NSET={boundary_set_name.upper()}"
    
    for line in lines:
        line_raw = line.strip()
        line_upper = line_raw.upper().replace(" ", "")
        
        if line_upper.startswith('*'):
            if target_header in line_upper:
                current_section = 'TARGET_NSET'
            else:
                current_section = None
            continue
        
        if current_section == 'TARGET_NSET':
            # 提取数字 ID，处理末尾逗号
            ids = [int(float(x)) for x in line_raw.rstrip(',').split(',')]
            boundary_node_ids.update(ids)

    print(f"提取到边界节点数量: {len(boundary_node_ids)}")

    # --- 3. 读取 velocity.txt 并准备插值 ---
    print("读取速度场数据并构建索引...")
    # 格式：x, y, vx, vy
    vel_data = np.loadtxt(velocity_file) 
    vel_coords = vel_data[:, :2] # x, y
    vel_values = vel_data[:, 2:] # vx, vy

    # 构建 k-d Tree
    tree = cKDTree(vel_coords)

    # 提取目标边界节点的坐标
    target_ids = sorted(list(boundary_node_ids))
    target_pts = np.array([nodes_dict[nid][:2] for nid in target_ids])

    # --- 4. 执行插值 (带外推逻辑的 IDW) ---
    print("执行速度场插值 (支持外推)...")
    # 寻找最近的 4 个点以进行加权平均，实现平滑外推
    distances, indices = tree.query(target_pts, k=4)

    # 处理距离为 0 的情况（点重合）并计算权重
    # 权重使用 1/d^2 (反距离平方加权)
    weights = 1.0 / np.maximum(distances, 1e-10)**2
    # 归一化权重
    weights_sum = np.sum(weights, axis=1)[:, np.newaxis]
    weights /= weights_sum

    # 计算插值后的 vx, vy (N x 2)
    # 分别对 vx 和 vy 进行加权
    interpolated_vel = np.sum(weights[:, :, np.newaxis] * vel_values[indices], axis=1)

    # --- 5. 导出 modelve.txt ---
    print("正在生成 modelve.txt ...")
    with open('modelve.txt', 'w', encoding='utf-8') as f:
        for i, nid in enumerate(target_ids):
            coords = nodes_dict.get(nid, [0.0, 0.0, 0.0])
            vx, vy = interpolated_vel[i]
            # 格式：节点编号, x, y, z, vx, vy
            f.write(f"{nid} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f} {vx:.6f} {vy:.6f}\n")

    print("速度边界处理完成。")


nodes, elements = parse_nodes_and_elements_fast('Job-1.inp')
set_results = parse_sets_from_inp('Job-1.inp')
filtered_el_dict = filter_element_nodes_by_z(nodes, elements, set_results)
node_neighbors = build_node_adjacency(filtered_el_dict, set_results['nodes'])
sorted_neighbors = sort_node_neighbors_custom(node_neighbors, nodes)
node_areas = calculate_neighbor_areas(sorted_neighbors, nodes)
node_attributes = identify_node_attributes(nodes, set_results['nodes'])
final_map, mapping_table = reindex_and_format_neighbors(set_results['nodes'], sorted_neighbors)


# 降雨插值程序（!!!precipitation.txt内容是x，y，pre；注意不是经纬度）
node_pre_results = {}
####node_pre_results = interpolate_precipitation(nodes, set_results['nodes'], 'precipitation.txt')

# 使用NC文件更新节点Z坐标并生成新INP文件(!!!terrain.nc文件应与模型边界1)
####update_nodes_z_with_nc('terrain.nc', nodes, set_results['nodes'], 'Job-1.inp', 'Job-1-new.inp')

# 插值速度边界条件（!!!velocity.txt内容是x，y，vx,vy；注意不是经纬度）
process_velocity_boundary('Job-1.inp', 'velocity-new.txt', nodes)

# 输出结果文件
output_model_files(final_map, node_areas, node_attributes, nodes, node_pre_results)

#执行读取
print(f"成功读取 {len(nodes)} 个节点")
for el_type, data in elements.items():
    print(f"单元类型 {el_type}: {len(data)} 个")

if set_results:
    print(f"节点集数量: {len(set_results['nodes'])}")
    print(f"单元集数量: {len(set_results['elements'])}")

# 打印结果示例
for el_type, el_map in filtered_el_dict.items():
    # 单元数：字典的键数量
    unit_count = len(el_map)
    # 去重后的节点数：把所有列表合并到集合
    unique_nodes = set().union(*el_map.values())
    node_count = len(unique_nodes)
    print(f"类型 {el_type} 过滤后剩余单元数: {unit_count}, 去重后剩余节点数: {node_count}")

#示例输出
####test_id=316
####print(f"节点 {test_id} 的连接节点有: {node_neighbors[test_id]}")
####print(f"节点 {test_id} 的原邻居: {node_neighbors[test_id]}")
####print(f"节点 {test_id} 的排序后邻居: {sorted_neighbors[test_id]}")
####print(f"节点 {test_id} 的计算面积为: {node_areas[test_id]}")

# 统计各属性数量
from collections import Counter
print(Counter(node_attributes.values()))

# 示例：查看前5个结果
####for nid in list(node_pre_results.keys())[:5]:
####    print(f"节点 {nid}: 插值 pre = {node_pre_results[nid]}")


####print(f"主节点{test_id} 对应的新编号邻居: {final_map[test_id]}")
