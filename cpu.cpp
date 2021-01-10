#include "cache.cpp"

class Cpu {

private:

    Cache * cache;
    int num_cache_sets;
    int block_size_doubles;
    int double_size_bytes;
    int instruction_count;

public:

    Cpu(Cache * cache, int num_cache_sets, int block_size_doubles, int double_size_bytes)
    {
        this->cache = cache;
        this->num_cache_sets = num_cache_sets;
        this->block_size_doubles = block_size_doubles;
        this->double_size_bytes = double_size_bytes;
        this->instruction_count = 0;
    }

    void reset_logging()
    {
        this->cache->reset_cache_logging();
        this->instruction_count = 0;
    }

    int get_read_hits()
    {
        return this->cache->get_read_hits();
    }

    int get_read_misses()
    {
        return this->cache->get_read_misses();
    }

    int get_write_hits()
    {
        return this->cache->get_write_hits();
    }

    int get_write_misses()
    {
        return this->cache->get_write_misses();
    }

    int get_instruction_count()
    {
        return this->instruction_count;
    }

    double loadDouble(Address * address)
    {
        this->instruction_count += 1;
        return this->cache->getDouble(address);
    }

    void storeDouble(Address * address, double value)
    {
        this->instruction_count += 1;
        this->cache->setDouble(address, value);
    }

    double addDouble(double value1, double value2)
    {
        this->instruction_count += 1;
        return value1 + value2;
    }

    double multDouble(double value1, double value2)
    {
        this->instruction_count += 1;
        return value1 * value2;
    }
};
