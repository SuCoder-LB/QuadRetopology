#include "decl.h"
#include "Pattern_all.h"

patchgen::VariableIndicators& patchgen::get_variable_indicators(int num_sides, int pattern_id) {
    typedef VariableIndicators& (*func_type)();
    static const func_type func_table_2[2] = {
        Pattern<2, 0>::get_variable_indicators,
        Pattern<2, 1>::get_variable_indicators
    };
    static const func_type func_table_3[4] = {
        Pattern<3, 0>::get_variable_indicators,
        Pattern<3, 1>::get_variable_indicators,
        Pattern<3, 2>::get_variable_indicators,
        Pattern<3, 3>::get_variable_indicators,
    };
    static const func_type func_table_4[5] = {
        Pattern<4, 0>::get_variable_indicators,
        Pattern<4, 1>::get_variable_indicators,
        Pattern<4, 2>::get_variable_indicators,
        Pattern<4, 3>::get_variable_indicators,
        Pattern<4, 4>::get_variable_indicators
    };
    static const func_type func_table_5[5] = {
        Pattern<5, 0>::get_variable_indicators,
        Pattern<5, 1>::get_variable_indicators,
        Pattern<5, 2>::get_variable_indicators,
        Pattern<5, 3>::get_variable_indicators,
        Pattern<5, 4>::get_variable_indicators,
    };
    static const func_type func_table_6[4] = {
        Pattern<6, 0>::get_variable_indicators,
        Pattern<6, 1>::get_variable_indicators,
        Pattern<6, 2>::get_variable_indicators,
        Pattern<6, 3>::get_variable_indicators,
    };
    static const func_type* func_table[7] = {
        nullptr,
        nullptr,
        func_table_2,
        func_table_3,
        func_table_4,
        func_table_5,
        func_table_6,
    };
    return func_table[num_sides][pattern_id]();
}
