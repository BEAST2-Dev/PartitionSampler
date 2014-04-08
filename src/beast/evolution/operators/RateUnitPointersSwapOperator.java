package beast.evolution.operators;

import beast.core.parameter.DPPointer;
import beast.evolution.alignment.GeneralUnitAlignment;
import beast.util.Randomizer;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 19/07/13
 * Time: 5:48 PM
 * To change this template use File | Settings | File Templates.
 */
public class RateUnitPointersSwapOperator extends RatePointersSwapOperator{
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
        int unit1 = sitesInCategoryIndex1[(int)(sitesInCategoryIndex1.length*Randomizer.nextDouble())];
        int unit2 = sitesInCategoryIndex2[(int)(sitesInCategoryIndex2.length*Randomizer.nextDouble())];

        int[] sitesInUnit1 = ((GeneralUnitAlignment)alignment).getSitesByUnit(unit1);
        int[] sitesInUnit2 = ((GeneralUnitAlignment)alignment).getSitesByUnit(unit2);

        boolean samePatterns = true;
        //System.out.println(sitesInUnit1.length +" "+ sitesInUnit2.length);
        if(sitesInUnit1.length == sitesInUnit2.length){

            for(int i = 0; i < sitesInUnit1.length; i++){
                if(alignment.getPatternIndex(sitesInUnit1[i]) != alignment.getPatternIndex(sitesInUnit2[i])){
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
        //System.out.println(samePatterns);

        DPPointer ratePointers = ratePointerInput.get(this);
        ratePointers.swapPointers(unit1, unit2);

        return 0.0;
    }
}
