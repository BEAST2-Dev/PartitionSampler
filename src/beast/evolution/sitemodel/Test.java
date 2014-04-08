package beast.evolution.sitemodel;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.QuietRealParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.SubstitutionModel;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 13/08/13
 * Time: 1:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class Test extends QuietSiteModel{
    public Input<QuietRealParameter> modelChoiceInput = new Input<QuietRealParameter>(
            "modelChoice",
            "And indicator that represents the current Gamma site model",
            Input.Validate.REQUIRED
    );

    public Input<Boolean> invPrLogitInput = new Input<Boolean>(
            "invPrLogit",
            "Whether invPr has logit transformation.",
            Input.Validate.REQUIRED
    );
    private boolean invPrLogit;
    private QuietRealParameter modelChoice;
    public void initAndValidate() throws Exception {
        this.modelChoice = modelChoiceInput.get();
        this.invPrLogit = invPrLogitInput.get();
        gammaCatCount = gammaCategoryCount.get();
        categoryCount = gammaCategoryCount.get();
        super.initAndValidate();

        addCondition(modelChoice);
    }

    private IntegerParameter indicator;

    public static final int SHAPE_INDEX = 0;
    public static final int INVAR_INDEX = 1;
    public static final int PRESENT = 1;
    public static final int ABSENT = 0;
    public static final int[][] INDICATORS = {
            {ABSENT, ABSENT,},
            {ABSENT, PRESENT},
            {PRESENT, ABSENT},
            {PRESENT, PRESENT}
        };



    public Test(){}
    public Test(
            SubstitutionModel substModel,
            QuietRealParameter muParameter,
            QuietRealParameter modelChoice,
            boolean invPr)throws Exception{

        this(substModel,
                muParameter,
                null,
                null,
                false,
                1,
                modelChoice,
                invPr);


    }

    public Test(SubstitutionModel substModel,
                             QuietRealParameter muParameter,
                             RealParameter shapeParameter,
                             RealParameter invarParameter,
                             boolean useBeast1StyleGamma,
                             int gammaCategoryCount,
                             QuietRealParameter modelChoice,
                             boolean invPrLogit) throws Exception{
        //m_bPropInvariantIsCategory = false;

        this.modelChoice = modelChoice;
        substitutionModel = (SubstitutionModel.Base)substModel;

        this.useBeast1StyleGamma = useBeast1StyleGamma;

        if (muParameter == null) {
            this.muParameter = new RealParameter("1.0");
        }else{
            this.muParameter = muParameter;
        }
        muParameter.setBounds(Math.max(muParameter.getLower(), 0.0), Math.min(muParameter.getUpper(), Double.POSITIVE_INFINITY));


        this.shapeParameter = shapeParameter;
        if (shapeParameter != null) {
            // The quantile calculator fails when the shape parameter goes much below
            // 1E-3 so we have put a hard lower bound on it. If this is not there then
            // the category rates can go to 0 and cause a -Inf likelihood (whilst this
            // is not a problem as the state will be rejected, it could mask other issues
            // and this seems the better approach.
            shapeParameter.setBounds(Math.max(shapeParameter.getLower(), 1.0E-3), Math.min(shapeParameter.getUpper(), 1.0E3));
        }

        this.invPrLogit = invPrLogit;
        if (invarParameter == null) {

            this.invarParameter = new RealParameter("0.0");
            this.invarParameter.setBounds(Math.max(0.0, this.invarParameter.getLower()), Math.min(1.0, this.invarParameter.getUpper()));

        }else{
            this.invarParameter = invarParameter;

        }



        if (/*invarParameter != null && */(this.invarParameter.getValue() < 0 || this.invarParameter.getValue() > 1) && !invPrLogit) {
            throw new Exception("proportion invariant should be between 0 and 1: "+this.invarParameter.getValue());
        }

        gammaCatCount = gammaCategoryCount;
        categoryCount = gammaCategoryCount;
        //System.out.println("gammaCategoryCount: "+gammaCategoryCount);
        //System.out.println("categoryCount: "+categoryCount);
        refresh();

        addCondition(muParameter);
        addCondition(invarParameter);
        addCondition(shapeParameter);
        addCondition(modelChoice);


    }

}
