package edu.albany.cs.base;

import org.apache.commons.lang3.ArrayUtils;

import java.util.Arrays;
import java.util.Comparator;

public class ArrayIndexComparator implements Comparator<Integer> {

    private final Double[] array;//we want to sort array by descending order
    public final Integer[] indexes;//keep indexes in descending order

    public ArrayIndexComparator(Double[] array) {
        this.array = array;
        this.indexes = createIndexArray();
    }

    public ArrayIndexComparator(double[] array) {
        this.array = ArrayUtils.toObject(array);
        this.indexes = createIndexArray();
    }

    public Integer[] createIndexArray() {
        Integer[] indexes = new Integer[array.length];
        for (int i = 0; i < this.array.length; i++) {
            indexes[i] = i;
        }
        return indexes;
    }

    @Override
    public int compare(Integer index1, Integer index2) {
        return array[index2].compareTo(array[index1]); // keep the descending order i1>i2> ...
    }

    public static void main(String args[]) {
        Double[] vectorRatioCB = new Double[]{0.1, 0.03, 0.2, 0.04, 0.5};
        ArrayIndexComparator arrayIndexComparator = new ArrayIndexComparator(vectorRatioCB);
        Integer[] indexes = arrayIndexComparator.indexes;
        Arrays.sort(indexes, arrayIndexComparator);
        System.out.println(Arrays.toString(indexes));
        for(int i:indexes){
            System.out.print(vectorRatioCB[i]+" ");
        }
    }
}