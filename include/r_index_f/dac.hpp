#ifndef _DAC_VECTOR_RIF_HH
#define _DAC_VECTOR_RIF_HH

#include <utility>
#include <iostream> 
#include <sdsl/int_vector.hpp>
#include <sdsl/rmq_support.hpp>
#include <vector>
#include <cmath>

using namespace sdsl;

class dac
{
public:
    typedef unsigned long int ulint;

    dac() {}

    dac(const std::vector<ulint> vec)
    {
        std::vector<bool> min_byte_rep = std::vector<bool>();
        std::vector<bool> bit_reps = std::vector<bool>();
        for(size_t i = 0; i<vec.size(); ++i)
        {
            ulint n = vec[i];
            int bytes = min_bytes(n);
            min_byte_rep.push_back(true);
            // Reduce by 1 since each needs at least one byte
            for(int b = 0; b < (bytes - 1); b++)
            {
                min_byte_rep.push_back(false);
            }

            ulint bits = bytes*8;
            // Gets the bit to be stored, from left to right
            while (bits>0) {
                bit_reps.push_back((n>>(bits-1))&1);
                --bits;
            }
        }

        byte_sizes = bv(min_byte_rep);
        select = bit_vector::select_1_type(&byte_sizes);
        values = bv(bit_reps);
    }

    //! []-operator
    ulint operator[](size_t i)
    {
        ulint pos_i = select(i + 1);
        ulint j = pos_i + 1;
        int curr_bytes = 1;
        while (j < byte_sizes.size() && (byte_sizes[j] != 1))
        {
            ++curr_bytes;
            ++j;
        }

        return get_int(pos_i*8, curr_bytes*8);
    }

    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "")
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        written_bytes += byte_sizes.serialize(out,v,"Byte_Sizes");
        written_bytes += values.serialize(out,v,"Values");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream &in)
    {
        byte_sizes.load(in);
        select = bit_vector::select_1_type(&byte_sizes);
        values.load(in);
    }

private:
    bit_vector byte_sizes;
    bit_vector::select_1_type select;
    bit_vector values;

    bit_vector bv(std::vector<bool> &b){

		if(b.size()==0) return bit_vector();

		bit_vector bv(b.size());

		for(uint64_t i=0;i<b.size();++i)
			bv[i] = b[i];

		return bv;
	}

    int min_bytes(ulint n)
    {
        int num_bytes = 0;
        do
        {
            n >>= 8;
            ++num_bytes;
        }
        while(n > 0);

        return num_bytes;
    }

    ulint get_int(size_t i, size_t len)
    {
        ulint num = 0;
        const ulint max_idx = len+i;
        while(i < max_idx)
        {
            num = (num << 1) | values[i++];
        }
        return num;
    }
};

#endif /* end of include guard: _DAC_VECTOR_RIF_HH */