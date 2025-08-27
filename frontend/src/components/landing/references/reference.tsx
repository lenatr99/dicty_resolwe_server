import React, { ReactElement, ReactNode } from 'react';
import { ReferenceContainer, ReferenceAuthors, ReferenceSource } from './references.styles';

type Props = {
    children: ReactNode;
    authors: string;
    source: string;
};

const Reference = ({ children, authors, source }: Props): ReactElement => (
    <ReferenceContainer>
        <p>{children}</p>
        <ReferenceAuthors>{authors}</ReferenceAuthors>
        <ReferenceSource>{source}</ReferenceSource>
    </ReferenceContainer>
);

export default Reference;
