from openmm.app import PDBReporter
import os



# 重写pdb文件保存逻辑
class DynamicPDBReporter(PDBReporter):
    def __init__(self, file_name, save_dir, interval):
        self.file_name = file_name
        self.save_dir = save_dir
        self.interval = interval
        super().__init__(os.path.join(save_dir, "placeholder.pdb"), interval)

    def report(self, simulation, state):
        current_step = simulation.context.getState().getStepCount()
        pdb_filename = os.path.join(self.save_dir, f"{self.file_name}_{current_step}.pdb")
        with open(pdb_filename, 'w') as f:
            self._out = f
            super().report(simulation, state)
        # self._out.close()
    
    def __del__(self):
        try:
            if not self._out.closed:  # 检查文件是否关闭
                self.writeFooter()
                self._out.close()
        except Exception as e:
            print(f"Error during cleanup: {e}")

        
if __name__ == "__main__":
    # 定义保存路径
    save_dir = "./output"
    os.makedirs(save_dir, exist_ok=True)
    # 添加 DynamicPDBReporter
    reporter = DynamicPDBReporter(save_dir, 1000)  # 每 1000 步保存一次
    simulation.reporters.append(reporter)