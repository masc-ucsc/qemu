/*
ESESC: Enhanced Super ESCalar simulator
Copyright (C) 2009 University of California, Santa Cruz.

Contributed by Jose Renau
Gabriel Southern

This file is part of ESESC.

ESESC is free software; you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation;
either version 2, or (at your option) any later version.

ESESC is    distributed in the  hope that  it will  be  useful, but  WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should  have received a copy of  the GNU General  Public License along with
ESESC; see the file COPYING.  If not, write to the  Free Software Foundation, 59
Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
using namespace std;

typedef uint32_t FlowID; // DInst.h

extern "C" void QEMUReader_goto_sleep(void *env);
extern "C" void QEMUReader_wakeup_from_sleep(void *env);
extern "C" void esesc_set_timing(uint32_t fid);

uint64_t last_pc = 0;
uint8_t ctrl_flag = 0;

uint8_t checkpc(uint64_t pc) {
  if (ctrl_flag == 1) {
    ctrl_flag = 0;
    if (pc != last_pc && last_pc != 0) {
      printf("error for ctrl jar pc: 0x%llx\n", (long long)last_pc);
      return 0;
    }
  } else if (pc == last_pc + 2 || pc == last_pc + 4 || last_pc == 0) {
    last_pc = pc;
    return 1;
  } else {
    printf("error pc increment: 0x%llx -> 0x%llx\n", (long long)last_pc,
           (long long)pc);
    last_pc = pc;
    return 0;
  }
  return 1;
}

string trans_reg(uint16_t src) {
  vector<string> register_list = {"",   "ra", "sp", "gp", "tp",
                                  "t0", "t1", "t2", "s0", "s1"};
  string res;
  if (src <= 9) {
    return (register_list[src]);
  } else if (src >= 10 && src <= 17) {
    res = "a" + to_string(src - 10);
    return (res);
  } else if (src >= 18 && src <= 27) {
    res = "s" + to_string(src - 16);
    return (res);
  } else if (src >= 28 && src <= 31) {
    res = "t" + to_string(src - 25);
    return (res);
  } else
    return (res);
}

extern "C" uint32_t QEMUReader_getFid(FlowID last_fid) { return 0; }

extern "C" void QEMUReader_finish_thread(uint32_t fid) {
  printf("QEMUReader_finishThread(%d)\n", fid);
}

extern "C" uint64_t QEMUReader_get_time() { return 0; }

extern "C" uint32_t QEMUReader_setnoStats() { return 0; }

extern "C" void QEMUReader_toggle_roi(uint32_t fid) {}

extern "C" uint64_t QEMUReader_queue_load(uint64_t pc, uint64_t addr,
                                          uint64_t data, uint16_t fid,
                                          uint16_t src1, uint16_t dest) {
  if (checkpc(pc))
    printf("ld   pc:0x%llx\t\tsrc1:%s\t\t\t\tdest:%s\t\taddr=0x%llx\n",
           (long long)pc, trans_reg(src1).c_str(), trans_reg(dest).c_str(),
           (long long)addr);
  return 0;
}
extern "C" uint64_t QEMUReader_queue_inst(uint64_t pc, uint64_t addr, int fid,
                                          int op, int src1, int src2, int dest,
                                          void *env) {
  // printf("%d pc=0x%llx addr=0x%llx op=%d src1=%d src2=%d dest=%d\n",fid,(long
  // long)pc,(long long)addr, op, src1, src2, dest);
  if (checkpc(pc))
    printf("alu  pc:0x%llx\top=%d\tsrc1:%s\t\tsrc2:%s\t\tdest:%s\n",
           (long long)pc, op, trans_reg(src1).c_str(), trans_reg(src2).c_str(),
           trans_reg(dest).c_str());
  return 0;
}

extern "C" uint64_t QEMUReader_queue_store(uint64_t pc, uint64_t addr,
                                           uint64_t data_new, uint64_t data_old,
                                           uint16_t fid, uint16_t src1,
                                           uint16_t src2, uint16_t dest) {
  if (checkpc(pc))
    printf("st   pc:0x%llx\t\tsrc1:%s\t\tsrc2:%s\t\t\t\taddr=0x%llx\n",
           (long long)pc, trans_reg(src1).c_str(), trans_reg(src2).c_str(),
           (long long)addr);
  return 0;
}

extern "C" uint64_t QEMUReader_queue_ctrl_data(uint64_t pc, uint64_t addr,
                                               uint64_t data1, uint64_t data2,
                                               uint16_t fid, uint16_t op,
                                               uint16_t src1, uint16_t src2,
                                               uint16_t dest) {
  // printf("pc=0x%llx addr=0x%llx last=0x%llx\n",(long long)pc,(long
  // long)addr,last_pc);
  if (checkpc(pc)) {
    printf("\nctrl pc:0x%llx\top=%d\tsrc1:%s\t\tsrc2:%s\t\tdest:%s\t\tjump to "
           "pc:0x%llx\n\n",
           (long long)pc, op, trans_reg(src1).c_str(), trans_reg(src2).c_str(),
           trans_reg(dest).c_str(), (long long)addr);
    last_pc = addr;
    ctrl_flag = 1;
  }
  return 0;
}

extern "C" void QEMUReader_finish(uint32_t fid) {
  printf("QEMUReader_finish(%d)\n", fid);
  exit(0);
}

extern "C" void QEMUReader_syscall(uint32_t num, uint64_t usecs, uint32_t fid) {
}

extern "C" FlowID QEMUReader_resumeThread(FlowID uid, FlowID last_fid) {
  static int conta = 0;
  printf("QEMUReader_resumeThread(%d,%d)\n", uid, last_fid);
  return conta++;
}

extern "C" void QEMUReader_pauseThread(FlowID fid) {
  printf("QEMUReader_pauseThread(%d)\n", fid);
}

extern "C" void QEMUReader_setFlowCmd(bool *flowStatus) {}

extern "C" uint32_t QEMUReader_cpu_start(uint32_t cpuid) {
  static int cpu_counter = 0;
  printf("QEMUReader_cpu_start(%d)\n", cpuid);
  return cpu_counter++;
}

extern "C" void QEMUReader_cpu_stop(uint32_t cpuid) {
  printf("QEMUReader_cpu_stop(%d)\n", cpuid);
}

extern "C" int qemuesesc_main(int argc, char **argv, char **envp);

void qemu_start(int argc, char **argv) {
  char **qargv = (char **)malloc(argc * sizeof(char **));

  qargv[0] = (char *)"qemu";
  for (int j = 1; j < argc; j++) {
    qargv[j] = strdup(argv[j]);
  }

  qemuesesc_main(argc, qargv, 0);
}

int main(int argc, char **argv) { qemu_start(argc, argv); }

// mkdir build
// PATH_TO_QEMU/configure --enable-esesc --target-list=riscv64-linux-user
//
// g++ PATH_TO_QEMU/qemuesesc_barebone.cpp -L ./riscv64-linux-user/
// ./riscv64-linux-user/libqemu.so
//
// ./a.out PATH_TO_ESESC/bins/riscv64/spec00_crafty.riscv64 <
// PATH_TO_ESESC/bins/inputs/crafty.in
//
