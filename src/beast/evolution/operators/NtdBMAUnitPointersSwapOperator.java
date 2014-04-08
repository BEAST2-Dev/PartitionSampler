package beast.evolution.operators;

import beast.core.parameter.DPPointer;
import beast.evolution.alignment.GeneralUnitAlignment;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
public class NtdBMAUnitPointersSwapOperator extends NtdBMAPointersSwapOperator{
    public void initAndValidate(){
        super.initAndValidate();
        if(!(alignment instanceof GeneralUnitAlignment)){
            throw new RuntimeException("Must use GeneralUnitAlignment.");
        }
    }
    public double proposal(){
        int categoryCount = dpVal.getCategoryCount();
        if(categoryCount == 1){
            return Double.NEGATIVE_INFINITY;
        }

        int categoryIndex1 = (int)(Randomizer.nextDouble()*categoryCount);
        int categoryIndex2 = categoryIndex1;
        while(categoryIndex1 == categoryIndex2){
            int something = (int)(Randomizer.nextDouble()*categoryCount);
            //System.out.println("Hi! "+categoryCount+" "+something);
            categoryIndex2 = something;
        }

        int[] sitesInCategoryIndex1 = dpVal.getClusterSites(categoryIndex1);
        int[] sitesInCategoryIndex2 = dpVal.getClusterSites(categoryIndex2);
        int unit1 = sitesInCategoryIndex1[Randomizer.nextInt(sitesInCategoryIndex1.length)];
        int unit2 = sitesInCategoryIndex2[Randomizer.nextInt(sitesInCategoryIndex2.length)];

        int[] sitesInUnit1 = ((GeneralUnitAlignment)alignment).getSitesByUnit(unit1);
        int[] sitesInUnit2 = ((GeneralUnitAlignment)alignment).getSitesByUnit(unit2);

        boolean samePatterns = true;
        if(sitesInUnit1.length == sitesInUnit2.length){
            for(int i = 0; i < sitesInUnit1.length; i++){
                if(sitesInUnit1[i] != sitesInUnit2[i]){
                    samePatterns = false;
                    break;
                }

            }
        }else{
            samePatterns = false;
        }

        if(samePatterns){
            return Double.NEGATIVE_INFINITY;
        }

        //System.out.println(categoryIndex1+" "+categoryIndex2);


        DPPointer paramPointers = paramPointersInput.get(this);
        paramPointers.swapPointers(unit1, unit2);
        DPPointer modelPointers = modelPointersInput.get(this);
        modelPointers.swapPointers(unit1, unit2);
        DPPointer freqPointers = freqPointersInput.get(this);
        freqPointers.swapPointers(unit1, unit2);

        return 0.0;
    }
}
