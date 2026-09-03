#!/bin/bash
OUT="rela.out"
# =====================================================
# 2. CONVERGENCIA SCF
# =====================================================

echo ""
echo "========================================"
echo "CONVERGENCIA SCF"
echo "========================================"

# Muestra las últimas iteraciones SCF

grep "scf:" $OUT | tail -10


# =====================================================
# 3. FUERZAS FINALES
# =====================================================

echo ""
echo "========================================"
echo "FUERZAS FINALES"
echo "========================================"

# Busca la fuerza máxima final

grep -i "Max" $OUT | grep -i "force"


# =====================================================
# 4. ESTRUCTURA RELAJADA
# =====================================================

echo ""
echo "========================================"
echo "ARCHIVOS DE ESTRUCTURA RELAJADA"
echo "========================================"

# Lista archivos importantes generados

ls *.XV 2>/dev/null
ls STRUCT_OUT 2>/dev/null
ls *.STRUCT_OUT 2>/dev/null


# =====================================================
# 5. MOSTRAR ÚLTIMAS LÍNEAS DEL OUTPUT
# =====================================================

echo ""
echo "========================================"
echo "FINAL DEL OUTPUT"
echo "========================================"

# Muestra las últimas líneas del cálculo

tail -30 $OUT


# =====================================================
# FIN
# =====================================================

echo ""
echo "Análisis completado"
