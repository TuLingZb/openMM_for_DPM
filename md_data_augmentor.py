"""
通过md生成相似的pdb结构,用于扩散模型训练
"""

"""
用于结构预测后批量relax structure

Done
"""

from openmm.app import *
# from openmm.app import PME, HBonds
from openmm import *
from openmm.unit import *
from sys import stdout
import os

from custom import DynamicPDBReporter

#-------------------------------------------准备模拟环境-------------------------------------------#

def em_pdb(pdb_dir,save_dir,pdb,molecule_force,constrain_atom_names=[]):
    """
    :param pdb_dir: pdb文件所在目录
    :param save_dir: 保存目录
    :param pdb: pdb文件名
    :param micromolecule_force: 自定义力场文件
    :param constrain_atom_names: 需要约束的原子名称列表
    """
    # 读取力场文件, 指定水模型文件
    # uo2_force = '/Users/zhangbai/Desktop/论文代码复现/NADIFF/openmm/force_field/uo2.xml'
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml', molecule_force)

    # 准备PDB文件
    # pdb = PDBFile(os.path.join(pdb_dir, f'{pdb_name}.pdb'))  # 替换为您的PDB文件

    # 检索pdb文件
    modeller = Modeller(pdb.topology, pdb.positions)

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

    # 创建 CompoundIntegrator 并添加积分器
    compound_integrator = CompoundIntegrator()
    compound_integrator.addIntegrator(emin_integrator)  # Index 0

    # 组装模拟环境
    simulation = Simulation(modeller.topology, system, compound_integrator)
    simulation.context.setPositions(modeller.positions)

    #--------------------------------constrain atoms--------------------------------#
    restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
    system.addForce(restraint)
    restraint.addGlobalParameter('k', 500.0*kilojoules_per_mole/nanometer) # 500-1000-2000
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')
    for atom in pdb.topology.atoms():
        if atom.name in constrain_atom_names:
            # system.setParticleMass(atom.index, 0*amu)
            restraint.addParticle(atom.index, pdb.positions[atom.index])

    #----------------------------------em----------------------------------#
    # 能量最小化阶段 (使用 index 0 的 VerletIntegrator)
    compound_integrator.setCurrentIntegrator(0)
    print("Performing energy minimization...")
    # 最大迭代次数maxIterations 
    # specify a tolerance for when the energy should be considered to have converged
    simulation.minimizeEnergy(tolerance=1*kilojoule/mole/nanometer)


    # 输出日志
    # simulation.reporters.append(PDBReporter(os.path.join(save_dir, "em_.pdb"), 1000)) # 1000 * 2fs = 2ps
    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
            potentialEnergy=True, temperature=True, volume=True))
    simulation.reporters.append(StateDataReporter(os.path.join(save_dir, "em_log.txt"), 1000, step=True,
            potentialEnergy=True, temperature=True, volume=True))

    #-----------------------------------保存----------------------------------#
    # 保存模拟结果
    print('Saving...')
    simulation.saveState(os.path.join(save_dir,'md_100ps.xml'))
    sim_state = simulation.context.getState(positions=True)
    # TODO 保存模拟状态方便续跑
    positions = sim_state.getPositions() # 可以保存到xml文件
    modeller = Modeller(simulation.topology, positions) #
    # 除去水分子
    modeller.deleteWater()
    # 去除水分子, 只保存结构
    with open(os.path.join(save_dir, 'md_fl_nowater.pdb'), 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print('Done')


if __name__ == '__main__':
    specify_pdb_name = "39E_M3_EM"
    simulate_results_dir = "/Users/zhangbai/Desktop/论文代码复现/NADIFF/openmm/simu_results/em"
    pdb_dir = "/Users/zhangbai/Desktop/论文代码复现/NADIFF/openmm/pdb_files"
    uo2_force = '/Users/zhangbai/Desktop/论文代码复现/NADIFF/openmm/force_field/uo2.xml'
    constrain_atom_names = ["N1","C2","N3","C4","C5","C6"]  # 指定需要约束的原子名称 核酸嘧啶环
    error_log = "/Users/zhangbai/Desktop/论文代码复现/NADIFF/openmm/simu_results/error_log.txt"
    if os.path.exists(error_log):
        # 文件存在则清楚文件内容
        with open(error_log, 'w') as f:
            f.write('')

    for pdb_file in os.listdir(pdb_dir):
        try:
            # if pdb_name == f'{specify_pdb_name}.pdb':
            if pdb_file.endswith(".pdb"):  # 检查文件扩展名是否是 .pdb
                # pdb_file_path = os.path.join(pdb_dir, pdb_file)
                # TODO 截取文件名称
                pdb_name = os.path.splitext(pdb_file)[0]  # 获取文件名（不含扩展名）
                # 创建PDB模拟文件夹
                save_dir = os.path.join(simulate_results_dir,pdb_name)
                if not os.path.exists(save_dir):
                    os.makedirs(save_dir)
                pdb = PDBFile(os.path.join(pdb_dir, f'{pdb_file}'))  # 替换为您的PDB文件
                em_pdb(pdb_dir, save_dir, pdb, uo2_force, constrain_atom_names)  # 调用操作函数
            
            with open(error_log, 'a') as f:
                f.write(f"{pdb_name}: Done\n")
        except Exception as e:
            print(f"Error processing file {pdb_name}: {e}")
            with open(error_log, 'a') as f:
                f.write(f"Error processing file {pdb_name}: {e}\n")
            continue