package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.*;
import beast.evolution.likelihood.DPSepTreeLikelihood;
import beast.evolution.likelihood.ExtendedSepTempWVTreeLikelihood;
import beast.evolution.sitemodel.DPNtdRateSepSiteModel;
import beast.math.distributions.ConditionalParametricDistribution;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

import java.util.HashMap;

/**
 * @author Chieh-Hsi Wu
 */
public class ComplexGammaBMASAMSPriorOperator extends Operator {

    public Input<DPPointer> ratesPointersInput = new Input<DPPointer>(
            "ratesPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );


    public Input<ParameterList> ratesListInput = new Input<ParameterList>(
            "ratesList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> alphaListInput = new Input<ParameterList>(
            "alphaList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> invPrListInput = new Input<ParameterList>(
            "invPrList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> siteModelListInput = new Input<ParameterList>(
            "siteModelList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );


    public Input<ExtendedSepTempWVTreeLikelihood> tempLikelihoodInput = new Input<ExtendedSepTempWVTreeLikelihood>(
            "tempLikelihood",
            "The temporary likelihood given the data at site i",
            Input.Validate.REQUIRED
    );

    public Input<DPValuable> dpValuableInput = new Input<DPValuable>(
            "dpVal",
            "reports the counts in each cluster",
            Input.Validate.REQUIRED
    );

    public Input<DPSepTreeLikelihood> dpTreeLikelihoodInput = new Input<DPSepTreeLikelihood>(
            "dpTreeLik",
            "Tree likelihood that handle DPP",
            Input.Validate.REQUIRED
    );


    public Input<ConditionalParametricDistribution> ratesBaseDistrInput = new Input<ConditionalParametricDistribution>(
            "ratesBaseDistr",
            "The base distribution of the overall subsitution rate of a site.",
            Input.Validate.REQUIRED
    );

    public Input<ConditionalParametricDistribution> alphaBaseDistrInput = new Input<ConditionalParametricDistribution>(
            "alphaBaseDistr",
            "The Dirichlet prior base distribution for alpha.",
            Input.Validate.REQUIRED
    );

    public Input<ConditionalParametricDistribution> invPrBaseDistrInput = new Input<ConditionalParametricDistribution>(
            "invPrBaseDistr",
            "The Dirichlet prior base distribution for invariant proportion.",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> siteModelBaseDistrInput = new Input<ParametricDistribution>(
            "siteModelBaseDistr",
            "The Dirichlet prior base distribution for site model"
    );

    public Input<Boolean> testCorrectInput = new Input<Boolean>(
            "testCorrect",
            "Whether to check the likelihood calculations are consistent.",
            false
    );



    private DPPointer ratesPointers;
    private ParameterList ratesList;
    private ParameterList alphaList;
    private ParameterList invPrList;
    private ParameterList siteModelList;
    private int pointerCount;
    private ConditionalParametricDistribution ratesBaseDistr;
    private ConditionalParametricDistribution alphaBaseDistr;
    private ConditionalParametricDistribution invPrBaseDistr;
    private ParametricDistribution siteModelBaseDistr;
    //private CompoundDirichletProcess dp;
    private ExtendedSepTempWVTreeLikelihood tempLikelihood;
    private DPSepTreeLikelihood dpTreeLikelihood;
    HashMap<Double,double[]> modelNetworkMap= new HashMap<Double,double[]>();
    private boolean testCorrect;
    public void initAndValidate(){
        testCorrect = testCorrectInput.get();
        alphaList = alphaListInput.get();
        invPrList = invPrListInput.get();
        ratesList = ratesListInput.get();
        siteModelList = siteModelListInput.get();

        ratesPointers = ratesPointersInput.get();

        pointerCount = ratesPointers.getDimension();



        //dp = dpInput.get();
        //List<ParametricDistribution> distrs = dp.getBaseDistributions();
        ratesBaseDistr = ratesBaseDistrInput.get();

        alphaBaseDistr = alphaBaseDistrInput.get();
        invPrBaseDistr = invPrBaseDistrInput.get();
        siteModelBaseDistr = siteModelBaseDistrInput.get();

        tempLikelihood = tempLikelihoodInput.get();
        dpTreeLikelihood = dpTreeLikelihoodInput.get();

        modelNetworkMap.put(1.0,new double[]{3.0});
        modelNetworkMap.put(2.0,new double[]{3.0});
        modelNetworkMap.put(3.0,new double[]{1.0,2.0,4.0});
        modelNetworkMap.put(4.0,new double[]{3.0,5.0});
        modelNetworkMap.put(5.0,new double[]{4.0});
        //System.out.println("is null? "+(modelNetworkMap.get(5.0) == null));

    }


    public double proposal(){
        double logq = 0.0;
        //Pick two indcies at random
        int index1 = Randomizer.nextInt(pointerCount);
        int index2 = index1;
        while(index2 == index1){
            index2 = Randomizer.nextInt(pointerCount);
        }


        int clusterIndex1 = ratesPointers.indexInList(index1,ratesList);
        int clusterIndex2 = ratesPointers.indexInList(index2,ratesList);

        //If the randomly draw sites are from the same cluster, perform a split-move.
        if(clusterIndex1 == clusterIndex2){

            int[] clusterSites = dpValuableInput.get().getClusterSites(clusterIndex1);
            //System.out.println("split: ");
            double temp = split(index1, index2,clusterIndex1,clusterSites);
            //System.out.println("split: "+temp);
            logq += temp;


        }else{
            //If the the two randomly drawn sites are not from the same cluster, perform a merge-move.

            int[] cluster1Sites = dpValuableInput.get().getClusterSites(clusterIndex1);
            int[] cluster2Sites = dpValuableInput.get().getClusterSites(clusterIndex2);

            //logq = merge(index1, index2,clusterIndex1,clusterIndex2,cluster1Sites,cluster2Sites);
            double temp = merge(
                    index1,
                    index2,
                    clusterIndex1,
                    clusterIndex2,
                    cluster1Sites,
                    cluster2Sites
            );

            //System.out.println("merge: "+temp);
            logq = temp;

        }
        return logq;
    }

    public QuietRealParameter getSample(ParametricDistribution distr, double upper, double lower) throws Exception{

        Double[][] sampleVals = distr.sample(1);

        QuietRealParameter sampleParameter = new QuietRealParameter(sampleVals[0]);
        sampleParameter.setUpper(upper);
        sampleParameter.setLower(lower);

        return sampleParameter;
    }

    public QuietRealParameter getSample(
            ConditionalParametricDistribution distr,
            RealParameter condition,
            double upper,
            double lower) throws Exception{

        Double[] sampleVal = distr.sample(1, (int)(double)condition.getValue())[0];
        QuietRealParameter sampleParameter = new QuietRealParameter(sampleVal);
        sampleParameter.setUpper(upper);
        sampleParameter.setLower(lower);
        return sampleParameter;
    }

    public double split(int index1, int index2, int clusterIndex, int[] initClusterSites){
        try{
            double logqSplit = 0.0;

            int clusterID = ratesList.getParameterIDNumber(clusterIndex);

            //Create a parameter by sampling from the prior
            QuietRealParameter newSiteModel = getSample(siteModelBaseDistr, siteModelList.getUpper(), siteModelList.getLower());
            QuietRealParameter newRates = getSample(ratesBaseDistr, newSiteModel, ratesList.getUpper(), ratesList.getLower());
            QuietRealParameter newAlpha = getSample(alphaBaseDistr, newSiteModel, alphaList.getUpper(), alphaList.getLower());
            QuietRealParameter newInvPr = getSample(invPrBaseDistr, newSiteModel, invPrList.getUpper(), invPrList.getLower());






            //Remove the index 1 and index 2 from the cluster
            int[] clusterSites = new int[initClusterSites.length -2];
            int k = 0;
            for(int i = 0 ; i < initClusterSites.length;i++){
                if(initClusterSites[i] != index1 && initClusterSites[i] != index2){
                    clusterSites[k++] = initClusterSites[i];
                }
            }
            //Form a new cluster with index 1
            ratesPointers.point(index1,newRates);

            //Shuffle the cluster_-{index_1,index_2} to obtain a random permutation
            Randomizer.shuffle(clusterSites);

            //Create the weight vector of site patterns according to the order of the shuffled index.
            /*int[] tempWeights = new int[tempLikelihood.m_data.get().getPatternCount()];
            int patIndex;
            for(int i = 0; i < clusterSites.length; i++){
                patIndex = tempLikelihood.m_data.get().getPatternIndex(clusterSites[i]);
                tempWeights[patIndex] = 1;
            }

            tempLikelihood.setPatternWeights(tempWeights);*/
            tempLikelihood.setupPatternWeightsFromSites(clusterSites);

            //Site log likelihoods in the order of the shuffled sites
            double[] logLik1 = tempLikelihood.calculateLogP(
                    newAlpha.getValue(),
                    newInvPr.getValue(),
                    newRates.getValue(),
                    newSiteModel.getValue(),
                    clusterSites
            );
            double[] logLik2 = new double[clusterSites.length];
            for(int i = 0; i < logLik2.length; i++){
                logLik2[i] = dpTreeLikelihood.getSiteLogLikelihood(
                        DPNtdRateSepSiteModel.RATES,
                        clusterID,
                        clusterSites[i]
                );
            }



            double[] lik1 = new double[logLik1.length];
            double[] lik2 = new double[logLik2.length];

            double maxLog;
            //scale it so it may be more accurate
            for(int i = 0; i < logLik1.length; i++){
                maxLog = Math.max(logLik1[i],logLik2[i]);
                //System.out.println(i+" "+logLik1[i]+" "+logLik2[i]);
                if(Math.exp(maxLog) < 1e-100){
                    if(maxLog == logLik1[i]){
                        lik1[i] = 1.0;
                        lik2[i] = Math.exp(logLik2[i] - maxLog);
                        //System.out.println(i+" "+lik1[i]+" "+lik2[i]);
                    }else{
                        lik1[i] = Math.exp(logLik1[i] - maxLog);
                        lik2[i] = 1.0;
                        //System.out.println(i+" "+lik1[i]+" "+lik2[i]);
                    }
                }else{

                    lik1[i] = Math.exp(logLik1[i]);
                    lik2[i] = Math.exp(logLik2[i]);

                }

            }

            /*for(int i = 0; i < logLik1.length; i++){
                lik1[i] = Math.exp(logLik1[i]);
                lik2[i] = Math.exp(logLik2[i]);
                //System.out.println(lik1[i]+" "+lik2[i]);
            }*/
            /*for(int i = 0; i < clusterSites.length;i++){
                System.out.println("clusterSites: "+clusterSites[i]);

            }
            System.out.println("index 1: "+index1+" index2: "+index2);*/

            int cluster1Count = 1;
            int cluster2Count = 1;

            int[] sitesInCluster1 = new int[initClusterSites.length];
            sitesInCluster1[0] = index1;

            //Assign members of the existing cluster (except for indice 1 and 2) randomly
            //to the existing and the new cluster
            double psi1, psi2, newClusterProb, draw;
            for(int i = 0;i < clusterSites.length; i++){



                psi1 = cluster1Count*lik1[i];
                psi2 = cluster2Count*lik2[i];
                newClusterProb = psi1/(psi1+psi2);
                draw = Randomizer.nextDouble();
                if(draw < newClusterProb){

                    //System.out.println("in new cluster: "+clusterSites[i]);
                    sitesInCluster1[cluster1Count] = clusterSites[i];
                    //paramPointers.point(clusterSites[i],newParam);
                    //modelPointers.point(clusterSites[i],newModel);
                    //freqsPointers.point(clusterSites[i],newFreqs);
                    //ratesPointers.point(clusterSites[i],newRates);
                    logqSplit += Math.log(newClusterProb);
                    cluster1Count++;
                }else{
                    logqSplit += Math.log(1.0-newClusterProb);
                    cluster2Count++;
                }

            }

            int newSiteModelInt = (int)(double)newSiteModel.getValue();
            logqSplit += ratesBaseDistr.calcLogP(newRates,newSiteModelInt)
                    + alphaBaseDistr.calcLogP(newAlpha,newSiteModelInt)
                    + invPrBaseDistr.calcLogP(newInvPr,newSiteModelInt)
                    + siteModelBaseDistr.calcLogP(newSiteModel)
            ;

            //Perform a split
            ratesList  = ratesListInput.get(this);
            alphaList  = alphaListInput.get(this);
            invPrList  = invPrListInput.get(this);

            siteModelList  = siteModelListInput.get(this);
            ratesPointers = ratesPointersInput.get(this);

            ratesList.splitParameter(clusterIndex,newRates);
            alphaList.splitParameter(clusterIndex,newAlpha);
            invPrList.splitParameter(clusterIndex,newInvPr);
            siteModelList.splitParameter(clusterIndex,newSiteModel);


            //Form a new cluster with index 1
                ratesPointers = ratesPointersInput.get(this);
                for(int i = 0 ; i < cluster1Count ;i++){
                    ratesPointers.point(sitesInCluster1[i],newRates);

                }
            return -logqSplit;


        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }

    public double merge(
            int index1,
            int index2,
            int clusterIndex1,
            int clusterIndex2,
            int[] cluster1Sites,
            int[] cluster2Sites){



        double logqMerge = 0.0;

        int cluster1ID = ratesList.getParameterIDNumber(clusterIndex1);
        int cluster2ID = ratesList.getParameterIDNumber(clusterIndex2);


        HashMap<Integer,Integer> siteMap = new HashMap<Integer, Integer>();

        //The value of the merged cluster will have that of cluster 2 before the merge.
        QuietRealParameter mergedRates = ratesList.getParameter(clusterIndex2);
        QuietRealParameter mergedAlpha = alphaList.getParameter(clusterIndex2);
        QuietRealParameter mergedInvPr = invPrList.getParameter(clusterIndex2);
        QuietRealParameter mergedSiteModel = siteModelList.getParameter(clusterIndex2);

        //Create a vector that combines the site indices of the two clusters
        int[] mergedClusterSites = new int[cluster1Sites.length+cluster2Sites.length-2];

        int k = 0;
        for(int i = 0; i < cluster1Sites.length;i++){

            if(cluster1Sites[i] != index1){
                // For all members that are not index 1,
                // record the cluster in which they have been before the merge,
                // and assign them to the combined vector.
                siteMap.put(cluster1Sites[i],clusterIndex1);
                mergedClusterSites[k++] = cluster1Sites[i];
            }
        }


        for(int i = 0; i < cluster2Sites.length;i++){
            //All members in cluster 2 remains in cluster2 so no new pointer assignments
            if(cluster2Sites[i] != index2){
                // For all members that are not index 2,
                // record the cluster in which they have been before the merge,
                // and assign them to the combined vector.
                siteMap.put(cluster2Sites[i],clusterIndex2);
                mergedClusterSites[k++] = cluster2Sites[i];
            }
        }




        try{

            // Create a weight vector of patterns to inform the temporary tree likelihood
            // which set of pattern likelihoods are to be computed.

            k = 0;
            int[] sCluster1Sites = new int[cluster1Sites.length - 1];
            for(int i = 0; i < cluster1Sites.length; i++){
                if(cluster1Sites[i] != index1){
                    sCluster1Sites[k++] = cluster1Sites[i];
                }
            }
            tempLikelihood.setupPatternWeightsFromSites(sCluster1Sites);
            double[] cluster1SitesCluster2ParamLogLik = tempLikelihood.calculateLogP(
                    mergedAlpha.getValue(),
                    mergedInvPr.getValue(),
                    mergedRates.getValue(),
                    mergedSiteModel.getValue(),
                    sCluster1Sites
            );


            k = 0;
            int[] sCluster2Sites = new int[cluster2Sites.length - 1];
            for(int i = 0; i < cluster2Sites.length; i++){
                if(cluster2Sites[i] != index2){
                    sCluster2Sites[k++] = cluster2Sites[i];
                }
            }
            tempLikelihood.setupPatternWeightsFromSites(sCluster2Sites);
            RealParameter removedRates = ratesList.getParameter(clusterIndex1);
            RealParameter removedAlpha = alphaList.getParameter(clusterIndex1);
            RealParameter removedInvPr = invPrList.getParameter(clusterIndex1);
            RealParameter removedSiteModel = siteModelList.getParameter(clusterIndex1);
            double[] cluster2SitesCluster1ParamLogLik = tempLikelihood.calculateLogP(
                    removedAlpha.getValue(),
                    removedInvPr.getValue(),
                    removedRates.getValue(),
                    removedSiteModel.getValue(),
                    sCluster2Sites
            );


            //System.out.println("populate logLik1:");
            double[] logLik1 = new double[mergedClusterSites.length];
            for(int i = 0; i < (cluster1Sites.length-1); i++){
                //System.out.println(clusterIndex1+" "+mergedClusterSites[i]);

                 logLik1[i] = dpTreeLikelihood.getSiteLogLikelihood(
                         DPNtdRateSepSiteModel.RATES,
                         cluster1ID,
                         mergedClusterSites[i]
                 );
            }
            System.arraycopy(cluster2SitesCluster1ParamLogLik,0,logLik1,cluster1Sites.length-1,cluster2SitesCluster1ParamLogLik.length);

            double[] logLik2 = new double[mergedClusterSites.length];
            System.arraycopy(cluster1SitesCluster2ParamLogLik,0,logLik2,0,cluster1SitesCluster2ParamLogLik.length);

            //System.out.println("populate logLik2:");
            for(int i = cluster1SitesCluster2ParamLogLik.length; i < logLik2.length; i++){
                //System.out.println(clusterIndex2+" "+mergedClusterSites[i-cluster1SitesCluster2ParamLogLik.length]);
                logLik2[i] = dpTreeLikelihood.getSiteLogLikelihood(
                        DPNtdRateSepSiteModel.RATES,
                        cluster2ID,
                        mergedClusterSites[i]
                );
            }



            double[] lik1 = new double[logLik1.length];
            double[] lik2 = new double[logLik2.length];

            //scale it so it may be more accuate
            double maxLog;
            //scale it so it may be more accurate
            for(int i = 0; i < logLik1.length; i++){
                maxLog = Math.max(logLik1[i],logLik2[i]);
                //System.out.println(i+" "+logLik1[i]+" "+logLik2[i]);
                if(Math.exp(maxLog) < 1e-100){
                    if(maxLog == logLik1[i]){
                        lik1[i] = 1.0;
                        lik2[i] = Math.exp(logLik2[i] - maxLog);
                        //System.out.println(i+" "+lik1[i]+" "+lik2[i]);
                    }else{
                        lik1[i] = Math.exp(logLik1[i] - maxLog);
                        lik2[i] = 1.0;
                        //System.out.println(i+" "+lik1[i]+" "+lik2[i]);
                    }
                }else{

                    lik1[i] = Math.exp(logLik1[i]);
                    lik2[i] = Math.exp(logLik2[i]);

                }

            }

            /*for(int i = 0; i < logLik1.length;i++){

                lik1[i] = Math.exp(logLik1[i]);
                lik2[i] = Math.exp(logLik2[i]);
                //System.out.println(lik1[i]+" "+lik2[i]);
            }*/

            //Create a set of indices for random permutation
            int[] shuffle = new int[mergedClusterSites.length];
            for(int i = 0; i < shuffle.length;i++){
                shuffle[i] = i;
            }
            Randomizer.shuffle(shuffle);

            int cluster1Count = 1;
            int cluster2Count = 1;
            int cluster;
            double psi1, psi2, cluster1Prob;
            for(int i = 0; i < mergedClusterSites.length;i++){

                cluster = siteMap.get(mergedClusterSites[shuffle[i]]);
                psi1 = cluster1Count*lik1[shuffle[i]];
                psi2 = cluster2Count*lik2[shuffle[i]];

                if(testCorrect)
                 testCorrectness(i,cluster,
                        clusterIndex1,clusterIndex2,shuffle, mergedClusterSites,
                         lik1,lik2);

                cluster1Prob = psi1/(psi1+psi2);
                //System.out.println(cluster1Prob);
                if(cluster == clusterIndex1){
                    logqMerge += Math.log(cluster1Prob);
                    cluster1Count++;

                }else if(cluster == clusterIndex2){
                    logqMerge += Math.log(1-cluster1Prob);
                    cluster2Count++;

                }else{
                    throw new RuntimeException("Something is wrong.");
                }


            }

            int removedSiteModelInt = (int)(double)removedSiteModel.getValue();
            logqMerge += alphaBaseDistr.calcLogP(removedAlpha,removedSiteModelInt)
                    + invPrBaseDistr.calcLogP(removedInvPr,removedSiteModelInt)
                    + ratesBaseDistr.calcLogP(removedRates,removedSiteModelInt)
                    + siteModelBaseDistr.calcLogP(removedSiteModel)
            ;

            if(logqMerge > Double.NEGATIVE_INFINITY){
                ratesList.mergeParameter(clusterIndex1,clusterIndex2);
                alphaList.mergeParameter(clusterIndex1,clusterIndex2);
                invPrList.mergeParameter(clusterIndex1,clusterIndex2);
                siteModelList.mergeParameter(clusterIndex1,clusterIndex2);
                for(int i = 0; i < cluster1Sites.length;i++){
                    //Point every member in cluster 1 to cluster 2
                    ratesPointers.point(cluster1Sites[i],mergedRates);


                }
            }
        }catch(Exception e){
            throw new RuntimeException(e);
        }

        return logqMerge;

    }






    public void testCorrectness(
            int i,
            int cluster,
            int clusterIndex1,
            int clusterIndex2,
            int[] shuffle,
            int[] mergedClusterSites,
            double[] lik1,
            double[] lik2) throws Exception{
        //System.out.println("Test correctness!");

        int[] tempWeights = new int[tempLikelihood.dataInput.get().getPatternCount()];
        tempWeights[tempLikelihood.dataInput.get().getPatternIndex(mergedClusterSites[shuffle[i]])] = 1;
        tempLikelihood.setPatternWeights(tempWeights);
        double temp1 = Math.exp(tempLikelihood.calculateLogP(
                alphaList.getParameter(clusterIndex1).getValue(),
                invPrList.getParameter(clusterIndex1).getValue(),
                ratesList.getParameter(clusterIndex1).getValue(),
                siteModelList.getParameter(clusterIndex1).getValue(),
                new int[]{mergedClusterSites[shuffle[i]]})[0]
        );
        double temp2 = Math.exp(tempLikelihood.calculateLogP(
                alphaList.getParameter(clusterIndex2).getValue(),
                invPrList.getParameter(clusterIndex2).getValue(),
                ratesList.getParameter(clusterIndex2).getValue(),
                siteModelList.getParameter(clusterIndex2).getValue(),
                new int[]{mergedClusterSites[shuffle[i]]})[0]
        );
        //System.out.println(temp1 +" "+ lik1[shuffle[i]] +" "+ temp2 +" "+ lik2[shuffle[i]]);
        if(temp1 != lik1[shuffle[i]] || temp2 != lik2[shuffle[i]]){
            System.out.println("temp1");
            System.out.println("shuffle_i: "+shuffle[i]);
            System.out.println("mergedClusterSites[shuffle]: "+mergedClusterSites[shuffle[i]]);
            System.out.println("cluster: "+cluster);
            System.out.println(+mergedClusterSites.length+" "+lik1.length);
            for(int j = 0; j < lik1.length;j++){
                System.out.println("merged lik1: "+mergedClusterSites[j]+" "+lik1[j]);
            }
            for(int j = 0; j < lik2.length;j++){
                System.out.println("merged lik2: "+mergedClusterSites[j]+" "+lik2[j]);
            }
            throw new RuntimeException(temp1+" "+lik1[shuffle[i]]+" "+temp2+" "+lik2[shuffle[i]]);

        }

    }
}
