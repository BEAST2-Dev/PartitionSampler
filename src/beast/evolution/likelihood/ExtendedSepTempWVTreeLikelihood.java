package beast.evolution.likelihood;

import beast.core.parameter.QuietRealParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.GeneralUnitAlignment;
import beast.evolution.sitemodel.DPNtdBMAGIBMASepSiteModel;
import beast.evolution.sitemodel.DummySiteModel;
import beast.evolution.sitemodel.QuietGammaSiteBMA;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 19/07/13
 * Time: 2:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class ExtendedSepTempWVTreeLikelihood extends SepTempWVTreeLikelihood{

    public double[] calculateLogP (
            RealParameter rateParameter,
            int[] siteIndex) throws Exception{
        GeneralUnitAlignment alignment = (GeneralUnitAlignment)dataInput.get();

        double[] logPs = new double[siteIndex.length];



            //Stores a list of subst models
            ArrayList<SwitchingNtdBMA> ntdBMAMap = new ArrayList<SwitchingNtdBMA>();
            //Stores the ID number of the substModel and pattern weights
            HashMap<Integer,int[]> clusterWeightsMap = new HashMap<Integer,int[]>();
            //Stores the ID number of the substModel and the sites with that model
            HashMap<Integer, ArrayList<Integer>> clusterSitesMap = new HashMap<Integer, ArrayList<Integer>>();

            int patternIndex;
            int modelID;

            for(int i = 0; i < siteIndex.length;i++){

                patternIndex = alignment.getPatternIndex(siteIndex[i]);
                //System.out.println((SwitchingNtdBMA)dpNtdRateSepSiteModel.getModel(siteIndex[i]));
                SwitchingNtdBMA ntdBMA = ((SwitchingNtdBMA)dpNtdRateSepSiteModel.getModel(alignment.getUnitBySite(siteIndex[i])));
                modelID = ntdBMA.getIDNumber();

                if(clusterWeightsMap.containsKey(modelID)){
                    //increment the weight of the pattern at site siteIndex[i] by 1
                    clusterWeightsMap.get(modelID)[patternIndex]++;
                    //add the site to the cluster
                    clusterSitesMap.get(modelID).add(i);

                }else{
                    //Create a new pattern weights vector
                    int[] patternWeights = new int[alignment.getPatternCount()];
                    //Increment the weight of the pattern at site siteIndex[i]
                    patternWeights[patternIndex]++;
                    //Add the new pattern weight vector to the HashMap
                    clusterWeightsMap.put(modelID,patternWeights);

                    //Add the new substModel to the array list
                    ntdBMAMap.add(ntdBMA);

                    //Create a new array list to store sites
                    ArrayList<Integer> clusterSites = new ArrayList<Integer>();
                    //Add site index to the corresponding list
                    clusterSites.add(i);
                    //Put the list in the HashMap
                    clusterSitesMap.put(modelID,clusterSites);
                }

            }

            int substModelCount = ntdBMAMap.size();
            int site;

            for(int i = 0; i < substModelCount;i++){

                m_siteModel.substModelInput.setValue(ntdBMAMap.get(i), m_siteModel);
                ((DummySiteModel)m_siteModel).getRateParameter().setValueQuietly(0,rateParameter.getValue());
                modelID = ntdBMAMap.get(i).getIDNumber();
                setPatternWeights(clusterWeightsMap.get(modelID));
                calculateLogP();

                ArrayList<Integer> arrayIndex = clusterSitesMap.get(modelID);
                for(Integer index:arrayIndex){
                    site = siteIndex[index];
                    logPs[index] = patternLogLikelihoods[dataInput.get().getPatternIndex(site)] ;
                }
            }


        return logPs;

    }

    public double[] calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int[] siteIndex){

        if(dpNtdRateSepSiteModel instanceof DPNtdBMAGIBMASepSiteModel){


            return calculateLogP2(
                    modelParameters,
                    modelCode,
                    freqs,
                    siteIndex
            );

        } else{
            return calculateLogP1(
                    modelParameters,
                    modelCode,
                    freqs,
                    siteIndex
            );
        }

    }

    public double[] calculateLogP1(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int[] siteIndex){

        GeneralUnitAlignment alignment = (GeneralUnitAlignment)dataInput.get();
        double[] logPs = new double[siteIndex.length];
        setModelParameterVals(modelParameters,modelCode,freqs);

        try{
            //Stores a list of rates
            ArrayList<QuietRealParameter> rates = new ArrayList<QuietRealParameter>();
            //Stores the ID number of the substModel and pattern weights
            HashMap<Integer,int[]> clusterWeightsMap = new HashMap<Integer,int[]>();
            //Stores the ID number of the substModel and the sites with that model
            HashMap<Integer, ArrayList<Integer>> clusterSitesMap = new HashMap<Integer, ArrayList<Integer>>();

            int patternIndex;
            int rateID;

            for(int i = 0; i < siteIndex.length;i++){

                patternIndex = alignment.getPatternIndex(siteIndex[i]);
                QuietRealParameter rate = dpNtdRateSepSiteModel.getRate(alignment.getUnitBySite(siteIndex[i]));
                rateID = rate.getIDNumber();

                if(clusterWeightsMap.containsKey(rateID)){
                    //increment the weight of the pattern at site siteIndex[i] by 1
                    clusterWeightsMap.get(rateID)[patternIndex]++;
                    //add the site to the cluster
                    clusterSitesMap.get(rateID).add(i);

                }else{
                    //Create a new pattern weights vector
                    int[] patternWeights = new int[alignment.getPatternCount()];
                    //Increment the weight of the pattern at site siteIndex[i]
                    patternWeights[patternIndex]++;
                    //Add the new pattern weight vector to the HashMap
                    clusterWeightsMap.put(rateID,patternWeights);

                    //Add the new substModel to the array list
                    rates.add(rate);

                    //Create a new array list to store sites
                    ArrayList<Integer> clusterSites = new ArrayList<Integer>();
                    //Add site index to the corresponding list
                    clusterSites.add(i);
                    //Put the list in the HashMap
                    clusterSitesMap.put(rateID,clusterSites);
                }

            }

            int rateCount = rates.size();
            int site;

            for(int i = 0; i < rateCount;i++){

                m_siteModel.substModelInput.setValue(substModel, m_siteModel);
                ((DummySiteModel)m_siteModel).getRateParameter().setValueQuietly(0,rates.get(i).getValue());
                rateID = rates.get(i).getIDNumber();
                setPatternWeights(clusterWeightsMap.get(rateID));
                calculateLogP();

                ArrayList<Integer> arrayIndex = clusterSitesMap.get(rateID);
                for(Integer index:arrayIndex){
                    site = siteIndex[index];
                    logPs[index] = patternLogLikelihoods[dataInput.get().getPatternIndex(site)] ;
                }
            }

        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return logPs;
    }

    public double[] calculateLogP2(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int[] siteIndex){
        /*System.out.println("modelParameters: "+modelParameters);
        System.out.println("modelCode: "+modelCode);
        System.out.println("freqs: "+freqs); */

        GeneralUnitAlignment alignment = (GeneralUnitAlignment)dataInput.get();
        double[] logPs = new double[siteIndex.length];
        setModelParameterVals(modelParameters,modelCode,freqs);

        try{
            //Stores a list of rates
            ArrayList<Double> rates = new ArrayList<Double>();
            ArrayList<Double> alpha = new ArrayList<Double>();
            ArrayList<Double> invPr = new ArrayList<Double>();
            ArrayList<Double> siteModelChoice = new ArrayList<Double>();
            ArrayList<Integer> rateIDs = new ArrayList<Integer>();

            //Stores the ID number of the substModel and pattern weights
            HashMap<Integer,int[]> clusterWeightsMap = new HashMap<Integer,int[]>();
            //Stores the ID number of the substModel and the sites with that model
            HashMap<Integer, ArrayList<Integer>> clusterSitesMap = new HashMap<Integer, ArrayList<Integer>>();

            int patternIndex;
            int rateID;


            DPNtdBMAGIBMASepSiteModel dpSiteModel = (DPNtdBMAGIBMASepSiteModel) dpNtdRateSepSiteModel;

            for(int i = 0; i < siteIndex.length;i++){

                patternIndex = alignment.getPatternIndex(siteIndex[i]);
                //QuietRealParameter rate = dpSiteModel.getRate(alignment.getUnitBySite(siteIndex[i]));
                rateID = dpSiteModel.getRateID(alignment.getUnitBySite(siteIndex[i]));

                if(clusterWeightsMap.containsKey(rateID)){
                    //increment the weight of the pattern at site siteIndex[i] by 1
                    clusterWeightsMap.get(rateID)[patternIndex]++;
                    //add the site to the cluster
                    clusterSitesMap.get(rateID).add(i);

                }else{
                    //Create a new pattern weights vector
                    int[] patternWeights = new int[alignment.getPatternCount()];
                    //Increment the weight of the pattern at site siteIndex[i]
                    patternWeights[patternIndex]++;
                    //Add the new pattern weight vector to the HashMap
                    clusterWeightsMap.put(rateID,patternWeights);

                    //Add the new substModel to the array list
                    double rateVal = dpSiteModel.getRateValue(alignment.getUnitBySite(siteIndex[i]));
                    double alphaVal = dpSiteModel.getAlphaValue(alignment.getUnitBySite(siteIndex[i]));
                    double invPrVal = dpSiteModel.getInvPrValue(alignment.getUnitBySite(siteIndex[i]));
                    double siteModelChoiceVal = dpSiteModel.getModelChoiceValue(alignment.getUnitBySite(siteIndex[i]));
                    rates.add(rateVal);
                    alpha.add(alphaVal);
                    invPr.add(invPrVal);
                    siteModelChoice.add(siteModelChoiceVal);
                    rateIDs.add(rateID);


                    //Create a new array list to store sites
                    ArrayList<Integer> clusterSites = new ArrayList<Integer>();
                    //Add site index to the corresponding list
                    clusterSites.add(i);
                    //Put the list in the HashMap
                    clusterSitesMap.put(rateID,clusterSites);
                }

            }

            int rateCount = rates.size();
            int site;

            for(int i = 0; i < rateCount;i++){

                m_siteModel.substModelInput.setValue(substModel, m_siteModel);
                ((QuietGammaSiteBMA)m_siteModel).setMuValueQuietly(rates.get(i));
                ((QuietGammaSiteBMA)m_siteModel).setShapeValueQuietly(alpha.get(i));
                ((QuietGammaSiteBMA)m_siteModel).setInvPrValueQuietly(invPr.get(i));
                ((QuietGammaSiteBMA)m_siteModel).setModelChoiceQuietly(siteModelChoice.get(i));
                ((QuietGammaSiteBMA)m_siteModel).setRatesKnown(false);

                /*System.out.println("rate: "+rates.get(i));
                System.out.println("alpha: "+alpha.get(i));
                System.out.println("invPr: "+invPr.get(i));
                System.out.println("siteModel: "+siteModelChoice.get(i));*/

                //rateID = rates.get(i).getIDNumber();
                setPatternWeights(clusterWeightsMap.get(rateIDs.get(i)));
                calculateLogP();

                //ArrayList<Integer> arrayIndex = clusterSitesMap.get(rateID);
                ArrayList<Integer> arrayIndex = clusterSitesMap.get(rateIDs.get(i));
                for(Integer index:arrayIndex){
                    site = siteIndex[index];
                    logPs[index] = patternLogLikelihoods[dataInput.get().getPatternIndex(site)] ;
                }
            }

        }catch(Exception e){
            throw new RuntimeException(e);

        }

        return logPs;
    }




    public double[] calculateLogP (
            double alpha,
            double invPr,
            double rate,
            double siteModelChoice,
            int[] siteIndex) throws Exception{

        GeneralUnitAlignment alignment = (GeneralUnitAlignment)dataInput.get();

        double[] logPs = new double[siteIndex.length];



            //Stores a list of subst models
            ArrayList<SwitchingNtdBMA> ntdBMAMap = new ArrayList<SwitchingNtdBMA>();
            //Stores the ID number of the substModel and pattern weights
            HashMap<Integer,int[]> clusterWeightsMap = new HashMap<Integer,int[]>();
            //Stores the ID number of the substModel and the sites with that model
            HashMap<Integer, ArrayList<Integer>> clusterSitesMap = new HashMap<Integer, ArrayList<Integer>>();

            int patternIndex;
            int modelID;

            for(int i = 0; i < siteIndex.length;i++){

                patternIndex = alignment.getPatternIndex(siteIndex[i]);
                //System.out.println((SwitchingNtdBMA)dpNtdRateSepSiteModel.getModel(siteIndex[i]));
                SwitchingNtdBMA ntdBMA = ((SwitchingNtdBMA)dpNtdRateSepSiteModel.getModel(alignment.getUnitBySite(siteIndex[i])));
                modelID = ntdBMA.getIDNumber();

                if(clusterWeightsMap.containsKey(modelID)){
                    //increment the weight of the pattern at site siteIndex[i] by 1
                    clusterWeightsMap.get(modelID)[patternIndex]++;
                    //add the site to the cluster
                    clusterSitesMap.get(modelID).add(i);

                }else{
                    //Create a new pattern weights vector
                    int[] patternWeights = new int[alignment.getPatternCount()];
                    //Increment the weight of the pattern at site siteIndex[i]
                    patternWeights[patternIndex]++;
                    //Add the new pattern weight vector to the HashMap
                    clusterWeightsMap.put(modelID,patternWeights);

                    //Add the new substModel to the array list
                    ntdBMAMap.add(ntdBMA);

                    //Create a new array list to store sites
                    ArrayList<Integer> clusterSites = new ArrayList<Integer>();
                    //Add site index to the corresponding list
                    clusterSites.add(i);
                    //Put the list in the HashMap
                    clusterSitesMap.put(modelID,clusterSites);
                }

            }

            int substModelCount = ntdBMAMap.size();
            int site;

            for(int i = 0; i < substModelCount;i++){

                m_siteModel.substModelInput.setValue(ntdBMAMap.get(i), m_siteModel);
                //((DummySiteModel)m_siteModel).getRateParameter().setValueQuietly(0,rateParameter.getValue());
                ((QuietGammaSiteBMA)m_siteModel).setMuValueQuietly(rate);
                ((QuietGammaSiteBMA)m_siteModel).setShapeValueQuietly(alpha);
                ((QuietGammaSiteBMA)m_siteModel).setInvPrValueQuietly(invPr);
                ((QuietGammaSiteBMA)m_siteModel).setModelChoiceQuietly(siteModelChoice);
                ((QuietGammaSiteBMA)m_siteModel).setRatesKnown(false);

                modelID = ntdBMAMap.get(i).getIDNumber();
                setPatternWeights(clusterWeightsMap.get(modelID));
                calculateLogP();

                ArrayList<Integer> arrayIndex = clusterSitesMap.get(modelID);
                for(Integer index:arrayIndex){
                    site = siteIndex[index];
                    logPs[index] = patternLogLikelihoods[dataInput.get().getPatternIndex(site)] ;
                }
            }


        return logPs;

    }
}
