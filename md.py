"""
用于常规MD模拟

需要添加离子
"""

from openmm.app import *
# from openmm.app import PME, HBonds
from openmm import *
from openmm.unit import *
from sys import stdout
import os

from custom import DynamicPDBReporter

#-------------------------------------------准备模拟环境-------------------------------------------#
# pdb_name = "5xm8_pb"
pdb_name = "39E_M3_EM"
simulate_results_dir = "/Users/zhangbai/Desktop/论文代码复现/NADIFF/openmm/simu_results"
pdb_dir = "/Users/zhangbai/Desktop/论文代码复现/NADIFF/openmm/pdb_files"
log_dir = "/Users/zhangbai/Desktop/论文代码复现/NADIFF/openmm/log_files"

save_dir = os.path.join(simulate_results_dir,pdb_name)
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    print(f"Folder '{save_dir}' created.")
else:
    print(f"Folder '{save_dir}' already exists.")


# 读取力场文件, 指定水模型文件
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml','/Users/zhangbai/Desktop/论文代码复现/NADIFF/openmm/force_field/uo2.xml')

# 准备PDB文件
pdb = PDBFile(os.path.join(pdb_dir, f'{pdb_name}.pdb'))  # 替换为您的PDB文件

# 检索pdb文件
modeller = Modeller(pdb.topology, pdb.positions)

# modeller.deleteWater() 

# 增加氢键 (需要力场来确定氢原子的位置, 选择指定 pH 下最常见的质子化状态)
modeller.addHydrogens(forcefield, pH=7.0)
# 如果是非标准元素加氢, 需要先自定义xml文件
# Modeller.loadHydrogenDefinitions('glycam-hydrogens.xml')

# 添加溶剂
modeller.addSolvent(
    forcefield, 
    boxSize=Vec3(5.0, 5.0, 5.0)*nanometers, # 指定盒子尺寸
    model='tip3p', # 指定水模型 'tip3p', 'spce', 'tip4pew', 'tip5p', 'swm4ndp'. 
    ionicStrength=0.15*molar, # 指定离子浓度
    positiveIon='Na+', # 指定阳离子 'Cs+' 、 'K+' 、 'Li+' 、 'Na+'和 'Rb+'
    negativeIon='Cl-',# 指定阴离子 'Cl-' 、 'Br-' 、 'F-'和 'I-'
)

#-------------------------------------------创建模拟系统-------------------------------------------#

# 创建系统
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.2*nanometers,
    constraints=HBonds,
    switchDistance=1.0*nanometers
)
# 设置积分器
 # 用于能量最小化
# emin_integrator = VerletIntegrator(0.001*picoseconds) 
emin_integrator = LangevinIntegrator(300*kelvin, 1.0/picosecond, 0.001*picoseconds) 
# 用于NVT阶段
nvt_integrator = LangevinIntegrator(300*kelvin, 1.0/picosecond, 0.002*picoseconds) 
# 用于NPT阶段 
npt_integrator = LangevinIntegrator(300*kelvin, 1.0/picosecond, 0.002*picoseconds)  

# 创建 CompoundIntegrator 并添加积分器
compound_integrator = CompoundIntegrator()
compound_integrator.addIntegrator(emin_integrator)  # Index 0
compound_integrator.addIntegrator(nvt_integrator)  # Index 1
compound_integrator.addIntegrator(npt_integrator)  # Index 2


# 组装模拟环境
simulation = Simulation(modeller.topology, system, compound_integrator)
simulation.context.setPositions(modeller.positions)

#-------------------------------------------能量最小化-------------------------------------------#
# 能量最小化阶段 (使用 index 0 的 VerletIntegrator)
compound_integrator.setCurrentIntegrator(0)
print("Performing energy minimization...")
# 最大迭代次数maxIterations 
# specify a tolerance for when the energy should be considered to have converged
simulation.minimizeEnergy(tolerance=1*kilojoule/mole/nanometer)
# simulation.minimizeEnergy()

# 输出日志
simulation.reporters.append(PDBReporter(os.path.join(save_dir, "em_.pdb"), 1000)) # 1000 * 2fs = 2ps
# simulation.reporters.append(DynamicPDBReporter("em", save_dir, 1000)) # 1000 * 2fs = 2ps
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True, volume=True))
simulation.reporters.append(StateDataReporter(os.path.join(save_dir, "em_log.txt"), 1000, step=True,
        potentialEnergy=True, temperature=True, volume=True))

#-------------------------------------------NVT 阶段-------------------------------------------#
# NVT 阶段 (使用 index 1 的 LangevinIntegrator)
compound_integrator.setCurrentIntegrator(1)
# simulation.context.setVelocitiesToTemperature(350*kelvin)
print("Running NVT equilibration...")
simulation.step(5000)  # 对应时间：5000*2fs = 10ps

# simulation.reporters.append(PDBReporter('nvt.pdb', 1000))
# simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
#         potentialEnergy=True, temperature=True, volume=True))
# simulation.reporters.append(StateDataReporter("nvt_log.txt", 100, step=True,
#         potentialEnergy=True, temperature=True, volume=True))

#-------------------------------------------NPT 阶段-------------------------------------------#
# NPT 阶段 (使用 index 2 的 LangevinIntegrator 并加入 Barostat)
compound_integrator.setCurrentIntegrator(2)
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
simulation.context.reinitialize(preserveState=True)
print("Running NPT equilibration...")
simulation.step(5000)  # 对应时间：5000*2fs = 10ps


#-------------------------------------------MD前准备-------------------------------------------#
# write trajectories
# simulation.reporters.append(PDBReporter('md.pdb', 1000)) # 1000 * 2fs = 2ps

# 保留第一帧
sim_state = simulation.context.getState(positions=True)
positions = sim_state.getPositions() # 可以保存到xml文件
PDBFile.writeFile(simulation.topology, positions, open(os.path.join(save_dir, 'md_f0.pdb'), 'w'))

simulation.reporters.append(XTCReporter(os.path.join(save_dir, 'md.xtc'), 1000))

# 日志输出项
simulation.reporters.append(StateDataReporter(
    os.path.join(save_dir, 'md_log.csv'), 1000, # 1000 * 2fs = 2ps
    step=True, # 当前时间步的索引)
    time=True,
    # progress=True, # 模拟已完成的百分比
    # remainingTime=True, # 完成模拟需要多长时间的估计
    totalSteps=100000 ,  # 用于估计模拟进度
    kineticEnergy=True, 
    potentialEnergy=True,
    totalEnergy=True,
    temperature=True,
    volume=True,
    # density=True, 
    # speed=True, 
    separator=' | ',  # csv 分隔符号
    )
 )
# 中间状态存档
#TODO
# # 设置checkpoint
# simulation.reporters.append(CheckpointReporter('checkpnt.chk', 5000))
# # 加载chkpnt
# simulation.loadCheckpoint('state.chk')

#-------------------------------------------MD 阶段-------------------------------------------#
print("Running MD equilibration...")
simulation.step(5000)  # 对应时间：50,00*2fs = 10ps

# 保存模拟结果
print('Saving...')
simulation.saveState(os.path.join(save_dir,'md_100ps.xml'))
sim_state = simulation.context.getState(positions=True)
# TODO 保存模拟状态方便续跑
positions = sim_state.getPositions() # 可以保存到xml文件
modeller = Modeller(simulation.topology, positions) #
# 除去水分子
modeller.deleteWater()
# 删除任意链、残基、原子和键
# modeller.delete(Atoms_obj_list) 
# 去除水分子, 只保存结构
with open(os.path.join(save_dir, 'md_fl_nowater.pdb'), 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print('Done')

