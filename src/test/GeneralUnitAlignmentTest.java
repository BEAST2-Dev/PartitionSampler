package test;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.GeneralUnitAlignment;
import beast.evolution.alignment.Sequence;
import junit.framework.TestCase;

/**
 * @author Chieh-Hsi Wu
 */
public class GeneralUnitAlignmentTest extends TestCase {



    public void test1() throws Exception{
        Sequence taxa1 = new Sequence("taxa1", "AGAAATATGTCTGAT");
        Sequence taxa2 = new Sequence("taxa2", "AGAAATATGTCTGAT");
        Sequence taxa3 = new Sequence("taxa3", "AGAAATATGTCTGAT");

        GeneralUnitAlignment data = new GeneralUnitAlignment();
        data.initByName(
                "sequence", taxa1,
                "sequence", taxa2,
                "sequence", taxa3,
                "dataType", "nucleotide",
                "unitDefinition", "site"
        );

        int[] siteToUnitMap = new int[]{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
        int[] unitToSitesMap = new int[]{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
        int siteCount = 15;
        for(int i = 0; i < siteCount; i++){

            assertEquals(data.getUnitBySite(i),siteToUnitMap[i]);
        }
    }

    public void test2() throws Exception{
        Sequence taxa1 = new Sequence("taxa1", "AGAAATATGTCTGAT");
        Sequence taxa2 = new Sequence("taxa2", "AGAAATATGTCTGAT");
        Sequence taxa3 = new Sequence("taxa3", "AGAAATATGTCTGAT");

        GeneralUnitAlignment data = new GeneralUnitAlignment();
        data.initByName(
                "sequence", taxa1,
                "sequence", taxa2,
                "sequence", taxa3,
                "dataType", "nucleotide",
                "unitDefinition", "codon"
        );

        int[] siteToUnitMap = new int[]{0,1,2,0,1,2,0,1,2,0,1,2,0,1,2};
        int[][] unitToSitesMap = new int[3][5];
        unitToSitesMap[0] = new int[]{0,3,6,9,12};
        unitToSitesMap[1] = new int[]{1,4,7,10,13};
        unitToSitesMap[2] = new int[]{2,5,8,11,14};

        int siteCount = 15;
        for(int i = 0; i < siteCount; i++){
            assertEquals(data.getUnitBySite(i),siteToUnitMap[i]);
        }

        int unitCount = 3;
        for(int i = 0; i < unitCount; i ++){
            int[] sites = data.getSitesByUnit(i);
            for(int j = 0; j < sites.length; j++){
                assertEquals(sites[j],unitToSitesMap[i][j]);
            }
        }
    }

    public void test3() throws Exception{
        Sequence taxa1 = new Sequence("taxa1", "AGAAATATGTC");
        Sequence taxa2 = new Sequence("taxa2", "AGAAATATGTC");
        Sequence taxa3 = new Sequence("taxa3", "AGAAATATGTC");

        GeneralUnitAlignment data = new GeneralUnitAlignment();
        data.initByName(
                "sequence", taxa1,
                "sequence", taxa2,
                "sequence", taxa3,
                "dataType", "nucleotide",
                "unitDefinition", "codon"
        );

        int[] siteToUnitMap = new int[]{0,1,2,0,1,2,0,1,2,0,1};
        int[][] unitToSitesMap = new int[3][5];
        unitToSitesMap[0] = new int[]{0,3,6,9};
        unitToSitesMap[1] = new int[]{1,4,7,10};
        unitToSitesMap[2] = new int[]{2,5,8};

        int siteCount = 11;
        for(int i = 0; i < siteCount; i++){
            assertEquals(data.getUnitBySite(i),siteToUnitMap[i]);
        }

        int unitCount = 3;
        for(int i = 0; i < unitCount; i ++){
            int[] sites = data.getSitesByUnit(i);
            for(int j = 0; j < sites.length; j++){
                assertEquals(sites[j],unitToSitesMap[i][j]);
            }
        }
    }


    public void test4() throws Exception{
        Sequence taxa1 = new Sequence("taxa1", "AGAAATATGTCTGAT");
        Sequence taxa2 = new Sequence("taxa2", "AGAAATATGTCTGAT");
        Sequence taxa3 = new Sequence("taxa3", "AGAAATATGTCTGAT");

        GeneralUnitAlignment data = new GeneralUnitAlignment();
        data.initByName(
                "sequence", taxa1,
                "sequence", taxa2,
                "sequence", taxa3,
                "dataType", "nucleotide",
                "unitDefinition", "codonByGene",
                "geneEndingSitePositions", "6 15"
        );

        int[] siteToUnitMap = new int[]{0,1,2,0,1,2,3,4,5,3,4,5,3,4,5};
        int[][] unitToSitesMap = new int[6][];
        unitToSitesMap[0] = new int[]{0,3};
        unitToSitesMap[1] = new int[]{1,4};
        unitToSitesMap[2] = new int[]{2,5};
        unitToSitesMap[3] = new int[]{6,9,12};
        unitToSitesMap[4] = new int[]{7,10,13};
        unitToSitesMap[5] = new int[]{8,11,14};

        int siteCount = 15;
        for(int i = 0; i < siteCount; i++){
            assertEquals(data.getUnitBySite(i),siteToUnitMap[i]);
        }

        int unitCount = 6;
        for(int i = 0; i < unitCount; i ++){
            int[] sites = data.getSitesByUnit(i);
            for(int j = 0; j < sites.length; j++){
                assertEquals(sites[j],unitToSitesMap[i][j]);
                //System.out.println(sites[j]+" "+unitToSitesMap[i][j]);
            }
        }
    }

    public void test5() throws Exception{
        Sequence taxa1 = new Sequence("taxa1", "AGAAATATGTCTGAT");
        Sequence taxa2 = new Sequence("taxa2", "AGAAATATGTCTGAT");
        Sequence taxa3 = new Sequence("taxa3", "AGAAATATGTCTGAT");

        GeneralUnitAlignment data = new GeneralUnitAlignment();
        data.initByName(
                "sequence", taxa1,
                "sequence", taxa2,
                "sequence", taxa3,
                "dataType", "nucleotide",
                "unitDefinition", "gene",
                "geneEndingSitePositions", "7 15"
        );

        int[] siteToUnitMap = new int[]{0,0,0,0,0,0,0,1,1,1,1,1,1,1,1};
        int[][] unitToSitesMap = new int[2][];
        unitToSitesMap[0] = new int[]{0,1,2,3,4,5,6};
        unitToSitesMap[1] = new int[]{7,8,9,10,11,12,13,14};

        int siteCount = 15;
        for(int i = 0; i < siteCount; i++){
            assertEquals(data.getUnitBySite(i),siteToUnitMap[i]);
        }

        int unitCount = 2;
        for(int i = 0; i < unitCount; i ++){
            int[] sites = data.getSitesByUnit(i);
            for(int j = 0; j < sites.length; j++){
                //System.out.println(sites.length);
                assertEquals(sites[j],unitToSitesMap[i][j]);
                //System.out.println(sites[j]+" "+unitToSitesMap[i][j]);
            }
        }
    }

    public void test6() throws Exception{
        Sequence taxa1 = new Sequence("taxa1", "AGAAATATGTCTGAT");
        Sequence taxa2 = new Sequence("taxa2", "AGAAATATGTCTGAT");
        Sequence taxa3 = new Sequence("taxa3", "AGAAATATGTCTGAT");

        GeneralUnitAlignment data = new GeneralUnitAlignment();
        data.initByName(
                "sequence", taxa1,
                "sequence", taxa2,
                "sequence", taxa3,
                "dataType", "nucleotide",
                "siteToUnitMap", "1 1 1 2 1 1 3 3 3 4 7 5 6 5 6"
        );

        int[] siteToUnitMap = new int[]{0,0,0,1,0,0,2,2,2,3,4,5,6,5,6};
        int[][] unitToSitesMap = new int[7][];
        unitToSitesMap[0] = new int[]{0,1,2,4,5};
        unitToSitesMap[1] = new int[]{3};
        unitToSitesMap[2] = new int[]{6,7,8};
        unitToSitesMap[3] = new int[]{9};
        unitToSitesMap[4] = new int[]{10};
        unitToSitesMap[5] = new int[]{11,13};
        unitToSitesMap[6] = new int[]{12,14};

        int siteCount = 15;
        for(int i = 0; i < siteCount; i++){
            assertEquals(data.getUnitBySite(i),siteToUnitMap[i]);
        }

        int unitCount = 2;
        for(int i = 0; i < unitCount; i ++){
            int[] sites = data.getSitesByUnit(i);
            for(int j = 0; j < sites.length; j++){
                //System.out.println(sites.length);
                assertEquals(sites[j],unitToSitesMap[i][j]);
                //System.out.println(sites[j]+" "+unitToSitesMap[i][j]);
            }
        }
    }



}
