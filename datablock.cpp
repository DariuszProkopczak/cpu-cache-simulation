class DataBlock {

public:

    int size;
    double * data;

    DataBlock(int size)
    {
        this->size = size;
        this->data = new double[this->size];
    }

};
