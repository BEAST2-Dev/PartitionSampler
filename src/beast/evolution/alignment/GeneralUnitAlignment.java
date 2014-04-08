package beast.evolution.alignment;

import beast.core.Input;
import beast.evolution.alignment.Alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * @author Chieh-Hsi Wu
 */
public class GeneralUnitAlignment extends Alignment {
    public Input<String> unitDefinitionInput = new Input<String>(
            "unitDefinition",
            "The definition of the assignment unit in the Dirichlet process mixture model."
    );

    public Input<String> siteToUnitMapInput = new Input<String>(
            "siteToUnitMap",
            "A string integers that maps from the site to unit.",
            Input.Validate.XOR,
            unitDefinitionInput
    );

    public Input<String> geneEndingSitePosInput = new Input<String>(
            "geneEndingSitePositions",
            "A string of integers that specifies the starting position of each gene in the alignment.",
            Input.Validate.OPTIONAL
    );

    public Input<Boolean> reorderInput = new Input<Boolean>(
            "reorder",
            "Whether order not to reorder the map",
            true
    );



    private int[] sitesToUnitsMap;
    private int[][] unitsToSitesMap;


    public void initAndValidate() throws Exception{
        super.initAndValidate();

        sitesToUnitsMap = new int[getSiteCount()];
        String unitDefinition = unitDefinitionInput.get();

        if(siteToUnitMapInput.get() != null){
            String[] mapStr = siteToUnitMapInput.get().trim().split("\\s+");
            HashMap<Integer,Integer> unitCounts = new HashMap<Integer, Integer>();
            ArrayList<Integer>[] arrayLists = (ArrayList<Integer>[])new ArrayList[getSiteCount()];

            int unit;
            int listCount = 0;
            for(int i = 0; i < mapStr.length; i++){
                unit = Integer.parseInt(mapStr[i]);


                if(!unitCounts.containsKey(unit)){
                    if(reorderInput.get()){
                        unitCounts.put(unit,listCount);
                    }else{
                        unitCounts.put(unit,unit);
                    }

                    System.err.println("Label "+unit+" is mapped to unit " +unitCounts.get(unit));
                    listCount++;
                }
                sitesToUnitsMap[i] = unitCounts.get(unit);
                if(arrayLists[unitCounts.get(unit)] == null)
                    arrayLists[unitCounts.get(unit)] = new ArrayList<Integer>();

                arrayLists[unitCounts.get(unit)].add(i);


            }

            unitsToSitesMap = new int[listCount][];
            for(int i = 0; i < unitsToSitesMap.length; i++){
                unitsToSitesMap[i] = new int[arrayLists[i].size()];
                for(int j = 0; j < unitsToSitesMap[i].length;j++){
                    unitsToSitesMap[i][j] = arrayLists[i].get(j);
                }
            }

        }else if(unitDefinition.equals("site")){


            sitesToUnitsMap = new int[getSiteCount()];
            unitsToSitesMap = new int[getSiteCount()][1];
            for(int i = 0; i < sitesToUnitsMap.length;i++){
                sitesToUnitsMap[i] = unitsToSitesMap[i][0] = i;
            }

        }else if(unitDefinition.equals("codon")){

            int[] geneEndingSitePos = new int[]{getSiteCount()};
            initiateCodonByGene(geneEndingSitePos);

        }else if(unitDefinition.equals("gene")){

            if(geneEndingSitePosInput.get() == null)
                throw new RuntimeException("If the unit is by gene/contigs, " +
                        "then the starting site positions of each gene/contigs.");

            String[] geneEndingSitePosStr = geneEndingSitePosInput.get().trim().split("\\s+");
            int[] geneEndingSitePos = new int[geneEndingSitePosStr.length];

            geneEndingSitePos[0] = Integer.parseInt(geneEndingSitePosStr[0]);
            for(int i = 1; i < geneEndingSitePos.length; i++){
                geneEndingSitePos[i] = Integer.parseInt(geneEndingSitePosStr[i]);
                if(geneEndingSitePos[i - 1] > geneEndingSitePos[i])
                    throw new RuntimeException("The ending positions need to be increasing order.");


            }

            if(geneEndingSitePos[geneEndingSitePos.length -1 ] != getSiteCount()){
                throw new RuntimeException("The last contig ending is not the end of the alignment.");
            }


            unitsToSitesMap = new int[geneEndingSitePos.length][];
            int prevEnd = 0;
            int k = 0;
            for(int i = 0; i < geneEndingSitePos.length; i++){
                unitsToSitesMap[i] = new int[geneEndingSitePos[i] - prevEnd];
                int length = geneEndingSitePos[i] - prevEnd;
                for(int j = 0; j < length; j++){
                    unitsToSitesMap[i][j] = k++;
                }
                prevEnd = geneEndingSitePos[i];
            }

            k = 0;
            prevEnd = 0;
            sitesToUnitsMap = new int[getSiteCount()];
            for(int i  = 0; i < geneEndingSitePos.length; i++){
                for(int j = prevEnd; j < geneEndingSitePos[i]; j++){
                    //System.out.println(i);
                    sitesToUnitsMap[k++] = i;
                }
                prevEnd = geneEndingSitePos[i];

            }


        }else if(unitDefinition.equals("codonByGene")){


            System.err.println("Assuming all genes start from the same codon position");

            if(geneEndingSitePosInput.get() == null)
                throw new RuntimeException("If the unit is by gene/contigs, " +
                        "then the starting site positions of each gene/contigs.");


            String[] geneEndingSitePosStr = geneEndingSitePosInput.get().trim().split("\\s+");
            int[] geneEndingSitePos = new int[geneEndingSitePosStr.length];
            geneEndingSitePos[0] = Integer.parseInt(geneEndingSitePosStr[0]);
            for(int i = 1; i < geneEndingSitePos.length; i++){
                geneEndingSitePos[i] = Integer.parseInt(geneEndingSitePosStr[i]);
                //System.out.println("Gene "+i+" ends at "+geneEndingSitePos[i]);
                if(geneEndingSitePos[i - 1] > geneEndingSitePos[i])
                    throw new RuntimeException("The ending positions need to be increasing order.");
            }

            initiateCodonByGene(geneEndingSitePos);



        }else{
            throw new RuntimeException("Please specify whether the assignment unit is site, genes, codon, codonByGene." +
                    "If it is none of the option above then please provide mapping of sites to units which should be a string of integers.");

        }

    }

    private void initiateCodonByGene(int[] geneEndingSitePos){
        for(int i =0; i < geneEndingSitePos.length; i++){
            System.out.println("gene " + i + " ends at " + geneEndingSitePos[i]);
        }


        if(geneEndingSitePos[geneEndingSitePos.length - 1] != getSiteCount()){
            throw new RuntimeException("The last contig ending is not the end of the alignment." +
                    "Site count is "+getSiteCount()+" but the end position of the last gene is "+geneEndingSitePos.length+".");
        }

        unitsToSitesMap = new int[geneEndingSitePos.length*3][];
        System.out.println(unitsToSitesMap.length);
        sitesToUnitsMap = new int[getSiteCount()];
        int k = 0;
        int prevEnd = 0;

        //Going through each gene segment
        for(int i = 0; i < geneEndingSitePos.length; i++){

            //Length of a geneSegment
            int length = geneEndingSitePos[i] - prevEnd;
            int dim = length/3;
            for(int j = 0; j < 3; j++){
                unitsToSitesMap[i*3+j] = new int[dim];
            }


            if(length  % 3 > 0){
                System.err.println("The number of sites in the gene is not a multiple of three.");
                dim += 1;
                for(int j = 0; j < length  % 3; j++){
                    unitsToSitesMap[i*3+j] = new int[dim];
                }
            }


            for(int j = 0; j < length; j++){
                //System.out.println(unitsToSitesMap[i*3 + j%3].length+" "+unitsToSitesMap.length);
                unitsToSitesMap[i*3 + j%3][j/3] = k;
                sitesToUnitsMap[k] = i*3 + j%3;
                k++;
            }

            prevEnd = geneEndingSitePos[i];

        }

    }

    public int getUnitBySite(int siteIndex){
        return sitesToUnitsMap[siteIndex];
    }

    public int getSiteCountGivenUnit(int unitIndex){
        return unitsToSitesMap[unitIndex].length;

    }

    public int[] getSitesByUnit(int unitIndex){
        return Arrays.copyOf(unitsToSitesMap[unitIndex], unitsToSitesMap[unitIndex].length);
    }



}
