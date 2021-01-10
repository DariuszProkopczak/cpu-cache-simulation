#include <cmath>

class Address {

private:

    int tag;
    int index;
    int offset;

public:

    int address;

    Address(int address, int num_cache_sets, int block_size, int bytes_per_value) {

        this->address = address;

        int num_bits_after_tag = (int) (log2(num_cache_sets) + log2(block_size) + log2(bytes_per_value));
        this->tag = address >> num_bits_after_tag;

        int num_bits_after_index = (int) (log2(block_size) + log2(bytes_per_value));
        int temp_address = address >> num_bits_after_index;
        this->index = temp_address % num_cache_sets;

        temp_address = address >> ((int) log2(bytes_per_value));
        int num_bits_offset = (int) log2(block_size);
        int offset_mask = (int) (pow(2, num_bits_offset) - 1);
        this->offset = temp_address & offset_mask;

    }

    int getTag()
    {
        return this->tag;
    }

    int getIndex()
    {
        return this->index;
    }

    int getOffset()
    {
        return this->offset;

    }
};
