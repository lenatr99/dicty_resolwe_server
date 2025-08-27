import React, { ReactElement, ReactNode } from 'react';
import { FeatureWrapper, FeatureDescription, FeatureImage } from './features.styles';

type Props = {
    children: ReactNode;
    name: string;
    imageSrc: string;
};

const Feature = ({ children, imageSrc, name }: Props): ReactElement => (
    <FeatureWrapper>
        <FeatureImage src={imageSrc} alt={name} />
        <FeatureDescription>{children}</FeatureDescription>
    </FeatureWrapper>
);

export default Feature;
