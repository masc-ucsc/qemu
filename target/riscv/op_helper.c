/*
 * RISC-V Emulation Helpers for QEMU.
 *
 * Copyright (c) 2016-2017 Sagar Karandikar, sagark@eecs.berkeley.edu
 * Copyright (c) 2017-2018 SiFive, Inc.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms and conditions of the GNU General Public License,
 * version 2 or later, as published by the Free Software Foundation.
 *
 * This program is distributed in the hope it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "qemu/osdep.h"
#include "qemu/log.h"
#include "cpu.h"
#include "qemu/main-loop.h"
#include "exec/exec-all.h"
#include "exec/helper-proto.h"

#ifdef CONFIG_ESESC
#include "esesc_qemu.h"
#include "InstOpcode.h"

#define AtomicAdd(ptr,val) __sync_fetch_and_add(ptr, val)
#define AtomicSub(ptr,val) __sync_fetch_and_sub(ptr, val)

volatile long long int icount = 0;
volatile long long int tcount = 0;

uint64_t esesc_mem_read(uint64_t addr);
uint64_t esesc_mem_read(uint64_t addr) {
  uint64_t buffer;

  CPUState *other_cs = first_cpu;

  CPU_FOREACH(other_cs) {
#ifdef CONFIG_USER_ONLY
    cpu_memory_rw_debug(other_cs, addr, &buffer, 8, 0);
#else
FIXME_NOT_DONE
    // FIXME: pass fid to mem_read and potentially call the system 
    // cpu_physical_memory_read(memaddr, myaddr, length);
#endif

    return buffer;
  }
  return 0;
}

void helper_esesc_dump(void);
void helper_esesc_dump(void) {
    CPUState *other_cs = first_cpu;

    CPU_FOREACH(other_cs) {
        printf("cpuid=%d halted=%d icount=%lld tcount=%lld\n",other_cs->fid, other_cs->halted,icount, tcount);
    }
}

void helper_esesc_load(CPURISCVState *env, uint64_t pc, uint64_t target, uint64_t data, uint64_t reg) {
  if (icount>0) {
    AtomicSub(&icount,1);
    return;
  }

  CPUState *cpu       = env_cpu(env);

  int src1 = reg & 0xFF;
  reg      = reg >> 8;
  reg      = reg >> 8;
  int dest = reg & 0xFF;
  AtomicAdd(&icount,QEMUReader_queue_load(pc, target, data, cpu->fid, src1, dest));
}

void helper_esesc_store(CPURISCVState *env, uint64_t pc, uint64_t target, uint64_t data_new, uint64_t data_old, uint64_t reg) {
  if (icount>0) {
    AtomicSub(&icount,1);
    return;
  }

  //printf("Hello pc:%llx data1:%lld data2:%lld\n",(long long)pc, (long long)data_new, (long long)data_old);

  CPUState *cpu       = env_cpu(env);

  int src1 = reg & 0xFF;
  reg      = reg >> 8;
  int src2 = reg & 0xFF;
  reg      = reg >> 8;
  int dest = reg & 0xFF;

  AtomicAdd(&icount,QEMUReader_queue_store(pc, target, data_new, data_old, cpu->fid, src1, src2, dest));
}

void helper_esesc_ctrl(CPURISCVState *env, uint64_t pc, uint64_t target, uint64_t op, uint64_t reg) {
  if (icount>0) {
    AtomicSub(&icount,1);
    return;
  }

  CPUState *cpu       = env_cpu(env);

  if (pc == target) {
    printf("jump to itself (terminate) pc:%llx\n",(unsigned long long)pc);
    QEMUReader_finish(cpu->fid);
    return;
  }

  int src1 = reg & 0xFF;
  reg      = reg >> 8;
  int src2 = reg & 0xFF;
  reg      = reg >> 8;
  int dest = reg & 0xFF;

  AtomicAdd(&icount,QEMUReader_queue_inst(pc, target, cpu->fid, op, src1, src2, dest));
}

void helper_esesc_ctrl_data(CPURISCVState *env, uint64_t pc, uint64_t target, uint64_t data1, uint64_t data2, uint64_t reg) {
  if (icount>0) {
    AtomicSub(&icount,1);
    return;
  }

  CPUState *cpu       = env_cpu(env);

  if (pc == target) {
    printf("jump to itself (terminate) pc:%llx\n",(long long)pc);
    QEMUReader_finish(cpu->fid);
    return;
  }

  int src1 = reg & 0xFF;
  reg      = reg >> 8;
  int src2 = reg & 0xFF;
  reg      = reg >> 8;
  int dest = reg & 0xFF;

#if 0
  if (pc==0x142cc)
    fprintf(stderr,"RAW: pc=%llx addr=%llx data1=%llx data2=%llx\n",(long long)pc,(long long)target,(long long)data1,(long long)data2);
#endif

  AtomicAdd(&icount,QEMUReader_queue_ctrl_data(pc, target, data1, data2, cpu->fid, iBALU_LBRANCH, src1, src2, dest));
}

void helper_esesc_alu(CPURISCVState *env, uint64_t pc, uint64_t op, uint64_t reg) {
  if (icount>0) {
    AtomicSub(&icount,1);
    return;
  }

  CPUState *cpu       = env_cpu(env);

  int src1 = reg & 0xFF;
  reg      = reg >> 8;
  int src2 = reg & 0xFF;
  reg      = reg >> 8;
  int dest = reg & 0xFF;

  AtomicAdd(&icount,QEMUReader_queue_inst(pc, 0, cpu->fid, op, src1, src2, dest));
}

#endif

/* Exceptions processing helpers */
void QEMU_NORETURN riscv_raise_exception(CPURISCVState *env,
                                          uint32_t exception, uintptr_t pc)
{
    CPUState *cs = env_cpu(env);
    qemu_log_mask(CPU_LOG_INT, "%s: %d\n", __func__, exception);
    cs->exception_index = exception;
    cpu_loop_exit_restore(cs, pc);
}

void helper_raise_exception(CPURISCVState *env, uint32_t exception)
{
    riscv_raise_exception(env, exception, 0);
}

target_ulong helper_csrrw(CPURISCVState *env, target_ulong src,
        target_ulong csr)
{
    target_ulong val = 0;
    if (riscv_csrrw(env, csr, &val, src, -1) < 0) {
        riscv_raise_exception(env, RISCV_EXCP_ILLEGAL_INST, GETPC());
    }
    return val;
}

target_ulong helper_csrrs(CPURISCVState *env, target_ulong src,
        target_ulong csr, target_ulong rs1_pass)
{
    target_ulong val = 0;
    if (riscv_csrrw(env, csr, &val, -1, rs1_pass ? src : 0) < 0) {
        riscv_raise_exception(env, RISCV_EXCP_ILLEGAL_INST, GETPC());
    }
    return val;
}

target_ulong helper_csrrc(CPURISCVState *env, target_ulong src,
        target_ulong csr, target_ulong rs1_pass)
{
    target_ulong val = 0;
    if (riscv_csrrw(env, csr, &val, 0, rs1_pass ? src : 0) < 0) {
        riscv_raise_exception(env, RISCV_EXCP_ILLEGAL_INST, GETPC());
    }
    return val;
}

#ifndef CONFIG_USER_ONLY

target_ulong helper_sret(CPURISCVState *env, target_ulong cpu_pc_deb)
{
    target_ulong prev_priv, prev_virt, mstatus;

    if (!(env->priv >= PRV_S)) {
        riscv_raise_exception(env, RISCV_EXCP_ILLEGAL_INST, GETPC());
    }

    target_ulong retpc = env->sepc;
    if (!riscv_has_ext(env, RVC) && (retpc & 0x3)) {
        riscv_raise_exception(env, RISCV_EXCP_INST_ADDR_MIS, GETPC());
    }

    if (get_field(env->mstatus, MSTATUS_TSR) && !(env->priv >= PRV_M)) {
        riscv_raise_exception(env, RISCV_EXCP_ILLEGAL_INST, GETPC());
    }

    mstatus = env->mstatus;

    if (riscv_has_ext(env, RVH) && !riscv_cpu_virt_enabled(env)) {
        /* We support Hypervisor extensions and virtulisation is disabled */
        target_ulong hstatus = env->hstatus;

        prev_priv = get_field(mstatus, MSTATUS_SPP);
        prev_virt = get_field(hstatus, HSTATUS_SPV);

        hstatus = set_field(hstatus, HSTATUS_SPV,
                                 get_field(hstatus, HSTATUS_SP2V));
        mstatus = set_field(mstatus, MSTATUS_SPP,
                            get_field(hstatus, HSTATUS_SP2P));
        hstatus = set_field(hstatus, HSTATUS_SP2V, 0);
        hstatus = set_field(hstatus, HSTATUS_SP2P, 0);
        mstatus = set_field(mstatus, SSTATUS_SIE,
                            get_field(mstatus, SSTATUS_SPIE));
        mstatus = set_field(mstatus, SSTATUS_SPIE, 1);

        env->mstatus = mstatus;
        env->hstatus = hstatus;

        if (prev_virt) {
            riscv_cpu_swap_hypervisor_regs(env);
        }

        riscv_cpu_set_virt_enabled(env, prev_virt);
    } else {
        prev_priv = get_field(mstatus, MSTATUS_SPP);

        mstatus = set_field(mstatus, MSTATUS_SIE,
                            get_field(mstatus, MSTATUS_SPIE));
        mstatus = set_field(mstatus, MSTATUS_SPIE, 1);
        mstatus = set_field(mstatus, MSTATUS_SPP, PRV_U);
        env->mstatus = mstatus;
    }

    riscv_cpu_set_mode(env, prev_priv);

    return retpc;
}

target_ulong helper_mret(CPURISCVState *env, target_ulong cpu_pc_deb)
{
    if (!(env->priv >= PRV_M)) {
        riscv_raise_exception(env, RISCV_EXCP_ILLEGAL_INST, GETPC());
    }

    target_ulong retpc = env->mepc;
    if (!riscv_has_ext(env, RVC) && (retpc & 0x3)) {
        riscv_raise_exception(env, RISCV_EXCP_INST_ADDR_MIS, GETPC());
    }

    target_ulong mstatus = env->mstatus;
    target_ulong prev_priv = get_field(mstatus, MSTATUS_MPP);
    target_ulong prev_virt = MSTATUS_MPV_ISSET(env);
    mstatus = set_field(mstatus, MSTATUS_MIE,
                        get_field(mstatus, MSTATUS_MPIE));
    mstatus = set_field(mstatus, MSTATUS_MPIE, 1);
    mstatus = set_field(mstatus, MSTATUS_MPP, PRV_U);
#ifdef TARGET_RISCV32
    env->mstatush = set_field(env->mstatush, MSTATUS_MPV, 0);
#else
    mstatus = set_field(mstatus, MSTATUS_MPV, 0);
#endif
    env->mstatus = mstatus;
    riscv_cpu_set_mode(env, prev_priv);

    if (riscv_has_ext(env, RVH)) {
        if (prev_virt) {
            riscv_cpu_swap_hypervisor_regs(env);
        }

        riscv_cpu_set_virt_enabled(env, prev_virt);
    }

    return retpc;
}

void helper_wfi(CPURISCVState *env)
{
    CPUState *cs = env_cpu(env);

    if ((env->priv == PRV_S &&
        get_field(env->mstatus, MSTATUS_TW)) ||
        riscv_cpu_virt_enabled(env)) {
        riscv_raise_exception(env, RISCV_EXCP_ILLEGAL_INST, GETPC());
    } else {
        cs->halted = 1;
        cs->exception_index = EXCP_HLT;
        cpu_loop_exit(cs);
    }
}

void helper_tlb_flush(CPURISCVState *env)
{
    CPUState *cs = env_cpu(env);
    if (!(env->priv >= PRV_S) ||
        (env->priv == PRV_S &&
         get_field(env->mstatus, MSTATUS_TVM))) {
        riscv_raise_exception(env, RISCV_EXCP_ILLEGAL_INST, GETPC());
    } else {
        tlb_flush(cs);
    }
}

void helper_hyp_tlb_flush(CPURISCVState *env)
{
    CPUState *cs = env_cpu(env);

    if (env->priv == PRV_M ||
        (env->priv == PRV_S && !riscv_cpu_virt_enabled(env))) {
        tlb_flush(cs);
        return;
    }

    riscv_raise_exception(env, RISCV_EXCP_ILLEGAL_INST, GETPC());
}

#endif /* !CONFIG_USER_ONLY */
