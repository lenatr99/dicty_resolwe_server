import { createContext } from 'react';
import { DictyRenderers } from 'redux/models/internal';

export const RendererContext = createContext<DictyRenderers>('canvas');
